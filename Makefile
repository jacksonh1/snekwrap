#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = snekwrap
PYTHON_VERSION = 3.10
PYTHON_INTERPRETER = python

#################################################################################
# COMMANDS                                                                      #
#################################################################################


## Install/update conda environment
.PHONY: requirements
requirements:
	mamba env update --name $(PROJECT_NAME) --file environment.yml


## Delete all compiled Python files ("*.py[co]", "__pycache__", etc.)
.PHONY: clean
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete


## create conda environment
.PHONY: create_environment
create_environment:
	mamba env create --name $(PROJECT_NAME) -f environment.yml

	@echo ">>> conda env created. Activate with:\nconda activate $(PROJECT_NAME)"


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################


## Format source code with black
.PHONY: format
format:
	black --config pyproject.toml snekwrap
#


# ## run script in data_processing to generate processed data
# .PHONY: data
# data:
# 	# $(PYTHON_INTERPRETER) data_processing/background_frequencies.py
# 	# $(PYTHON_INTERPRETER) data_processing/make_count_matrices.py


#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys; \
lines = '\n'.join([line for line in sys.stdin]); \
matches = re.findall(r'\n## (.*)\n[\s\S]+?\n([a-zA-Z_-]+):', lines); \
print('Available rules:\n'); \
print('\n'.join(['{:25}{}'.format(*reversed(match)) for match in matches]))
endef
export PRINT_HELP_PYSCRIPT

help:
	@$(PYTHON_INTERPRETER) -c "${PRINT_HELP_PYSCRIPT}" < $(MAKEFILE_LIST)


## Lint using flake8 and black (use `make format` to do formatting)
# .PHONY: lint
# lint:
# 	flake8 lir_proteome_screen_pssm
# 	isort --check --diff --profile black lir_proteome_screen_pssm
# 	black --check --config pyproject.toml lir_proteome_screen_pssm


## Format source code with black
# .PHONY: format
# format:
# 	black --config pyproject.toml lir_proteome_screen_pssm
#

## run scripts in data_processing to generate processed data
# .PHONY: data
# data: requirements
# 	$(PYTHON_INTERPRETER) data_processing/background_frequencies.py
# 	# $(PYTHON_INTERPRETER) data_processing/make_count_matrices.py
# 	$(PYTHON_INTERPRETER) data_processing/lir_central_test_set.py
# 	$(PYTHON_INTERPRETER) data_processing/lir_central_augmented_test_set.py
# 	$(PYTHON_INTERPRETER) data_processing/process_tables.py
