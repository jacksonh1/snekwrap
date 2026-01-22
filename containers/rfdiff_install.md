```bash
apptainer build rf_diffusion.sif ./rf_diff.apptainer_w_weights.def
singularity run --nv /home/jhalpin/orcd/pool/09-fragfold/tools/snekwrap/containers/rf_diffusion.sif --help
```









# ```
# mamba create -n SE3nv python=3.9
# conda activate SE3nv
# mamba install pytorch=2.0.1 torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia
# pip install dgl -f https://data.dgl.ai/wheels/cu121/repo.html
# mamba install hydra-core omegaconf icecream pyrsistent
# cd env/SE3Transformer
# pip install --no-cache-dir -r requirements.txt
# python setup.py install
# cd ../..
# pip install -e .
# ```
# 
# 
# The above was still fucked up even though it installed.
# 
# 
# ```
# docker build -f docker/Dockerfile -t rfdiffusion .
# docker save cfa17d1866d1 -o rfdimage.tar # ran docker images to find tag
# sudo singularity build rfd.sif docker-archive://rfdimage.tar
# ```
# 
# 
# That didn't work either because of permission issues. So I added a chmod line to the Dockerfile.fixed and built that one instead.
# 
# ```
# cd /mnt/shared2/jch/09-fragfold/tools/_third_party/RFdiffusion && docker build -f docker/Dockerfile.fixed -t rfdiffusion:fixed .
# cd /mnt/shared2/jch/09-fragfold/tools/_third_party/RFdiffusion && docker save rfdiffusion:fixed -o /mnt/shared2/jch/09-fragfold/tools/_third_party/RFdiffusion/rfdiffusion_fixed.tar
# cd /mnt/shared2/jch/09-fragfold/tools/_third_party/RFdiffusion && singularity build rfd_fixed.sif docker-archive://rfdiffusion_fixed.tar
# ```

