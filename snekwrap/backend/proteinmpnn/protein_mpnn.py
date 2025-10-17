"""
adapted from:
    https://colab.research.google.com/drive/1OU0KJegdWYZmEzDxPm73w14WkLSift6t
"""


import re
import numpy as np
import copy
import torch
import snekwrap.config as config
from pathlib import Path
from typing import Literal
from snekwrap.backend.proteinmpnn.protein_mpnn_utils import (
    tied_featurize, parse_PDB, StructureDatasetPDB, ProteinMPNN, _scores, _S_to_seq
)

class ProteinMPNNSample:
    def __init__(
        self,
        sequence: str,
        score: float,
        sample_number: int,
        input_params: dict,
    ):
        self.sequence = sequence
        self.score = score
        self.sample_number = sample_number
        self.input_params = input_params

def run_protein_mpnn(
    pdb_path: str | Path,
    designed_chains: list[str] | None = None,
    fixed_chains: list[str] | None = None,
    fixed_position_chain: str | None = None,
    fixed_positions: list[int] | None = None,
    num_seqs: int = 8,
    sampling_temp: str = "0.1",
    model_weights_folder: Literal["vanilla_model_weights", "soluble_model_weights"] = "soluble_model_weights",  # Only allow these two values
    backbone_noise: float = 0.20,
    model_name: str = "v_48_020",
    proteinmpnn_repo: str | Path = config.PROTEINMPNN_REPO,
) -> list[ProteinMPNNSample]:
    proteinmpnn_repo = Path(proteinmpnn_repo)
    if designed_chains is None:
        designed_chains = ["A"]
    device = torch.device("cuda:0" if (torch.cuda.is_available()) else "cpu")
    PATH_TO_MODEL_WEIGHTS = proteinmpnn_repo / model_weights_folder
    _hidden_dim = 128
    _num_layers = 3
    model_folder_path = PATH_TO_MODEL_WEIGHTS
    checkpoint_path = model_folder_path / f"{model_name}.pt"
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model = ProteinMPNN(
        num_letters=21,
        node_features=_hidden_dim,
        edge_features=_hidden_dim,
        hidden_dim=_hidden_dim,
        num_encoder_layers=_num_layers,
        num_decoder_layers=_num_layers,
        augment_eps=backbone_noise,
        k_neighbors=checkpoint["num_edges"],
    )
    model.to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    pdb_name = Path(pdb_path).stem

    # Accept lists directly for designed_chain, fixed_chain, and fixed_positions
    designed_chain_list = designed_chains if designed_chains is not None else []
    fixed_chain_list = fixed_chains if fixed_chains is not None else []
    chain_list = list(set(designed_chain_list + fixed_chain_list))

    # Check that fixed_position_chain is in designed_chain if defined
    if fixed_position_chain is not None and fixed_position_chain not in designed_chain_list:
        raise ValueError(f"fixed_position_chain '{fixed_position_chain}' must be present in 'designed_chain'. If you want to fix all positions in a chain, add it to 'fixed_chain' instead.")
    if fixed_position_chain is None:
        fixed_position_chain = ""

    if fixed_positions is not None and len(fixed_positions) > 0:
        fixed_positions_dict = {}
        fixed_positions_dict[pdb_name] = {}
        for c in chain_list:
            if c == fixed_position_chain:
                fixed_positions_dict[pdb_name][c] = [int(item) for item in fixed_positions]
            else:
                fixed_positions_dict[pdb_name][c] = []
    else:
        fixed_positions_dict = None

    batch_size = 1
    max_length = 20000
    num_seq_per_target = num_seqs
    NUM_BATCHES = num_seq_per_target // batch_size
    BATCH_COPIES = batch_size
    temperatures = [float(item) for item in str(sampling_temp).split()]
    omit_AAs = "X"
    alphabet = "ACDEFGHIKLMNPQRSTVWYX"
    omit_AAs_np = np.array([AA in omit_AAs for AA in alphabet]).astype(np.float32)
    bias_AAs_np = np.zeros(len(alphabet))
    chain_id_dict = {}
    chain_id_dict[pdb_name] = (designed_chain_list, fixed_chain_list)
    pssm_dict = None
    omit_AA_dict = None
    bias_AA_dict = None
    tied_positions_dict = None
    bias_by_res_dict = None
    pdb_path = str(pdb_path)
    pdb_dict_list = parse_PDB(pdb_path, input_chain_list=chain_list)
    dataset_valid = StructureDatasetPDB(pdb_dict_list, truncate=None, max_length=max_length)

    results = []
    with torch.no_grad():
        for ix, protein in enumerate(dataset_valid):
            batch_clones = [copy.deepcopy(protein) for _ in range(BATCH_COPIES)]
            (
                X, S, mask, lengths, chain_M, chain_encoding_all, chain_list_list,
                visible_list_list, masked_list_list, masked_chain_length_list_list,
                chain_M_pos, omit_AA_mask, residue_idx, dihedral_mask,
                tied_pos_list_of_lists_list, pssm_coef, pssm_bias, pssm_log_odds_all,
                bias_by_res_all, tied_beta,
            ) = tied_featurize(
                batch_clones, device, chain_id_dict, fixed_positions_dict,
                omit_AA_dict, tied_positions_dict, pssm_dict, bias_by_res_dict,
            )
            pssm_log_odds_mask = (pssm_log_odds_all > 0.0).float()
            name_ = batch_clones[0]["name"]
            randn_1 = torch.randn(chain_M.shape, device=X.device)
            log_probs = model(
                X, S, mask, chain_M * chain_M_pos, residue_idx, chain_encoding_all, randn_1
            )
            mask_for_loss = mask * chain_M * chain_M_pos
            scores = _scores(S, log_probs, mask_for_loss)
            native_score = scores.cpu().data.numpy()
            for temp in temperatures:
                for j in range(NUM_BATCHES):
                    randn_2 = torch.randn(chain_M.shape, device=X.device)
                    sample_dict = model.sample(
                        X,
                        randn_2,
                        S,
                        chain_M,
                        chain_encoding_all,
                        residue_idx,
                        mask=mask,
                        temperature=temp,
                        omit_AAs_np=omit_AAs_np,
                        bias_AAs_np=bias_AAs_np,
                        chain_M_pos=chain_M_pos,
                        omit_AA_mask=omit_AA_mask,
                        pssm_coef=pssm_coef,
                        pssm_bias=pssm_bias,
                        pssm_multi=0.0,
                        pssm_log_odds_flag=False,
                        pssm_log_odds_mask=pssm_log_odds_mask,
                        pssm_bias_flag=False,
                        bias_by_res=bias_by_res_all,
                    )
                    S_sample = sample_dict["S"]
                    log_probs = model(
                        X,
                        S_sample,
                        mask,
                        chain_M * chain_M_pos,
                        residue_idx,
                        chain_encoding_all,
                        randn_2,
                        use_input_decoding_order=True,
                        decoding_order=sample_dict["decoding_order"],
                    )
                    mask_for_loss = mask * chain_M * chain_M_pos
                    scores = _scores(S_sample, log_probs, mask_for_loss)
                    scores = scores.cpu().data.numpy()
                    for b_ix in range(BATCH_COPIES):
                        masked_chain_length_list = masked_chain_length_list_list[b_ix]
                        masked_list = masked_list_list[b_ix]
                        seq = _S_to_seq(S_sample[b_ix], chain_M[b_ix])
                        score = scores[b_ix]
                        input_params = dict(
                            PDB_PATH=pdb_path,
                            designed_chain=designed_chains,
                            fixed_chain=fixed_chains,
                            fixed_position_chain=fixed_position_chain,
                            fixed_positions=fixed_positions,
                            num_seqs=num_seqs,
                            sampling_temp=sampling_temp,
                            MODEL_WEIGHTS_FOLDER=model_weights_folder,
                            BACKBONE_NOISE=backbone_noise,
                            MODEL_NAME=model_name,
                        )
                        results.append(
                            ProteinMPNNSample(
                                sequence=seq,
                                score=score,
                                sample_number=j,
                                input_params=input_params,
                            )
                        )
    return results
