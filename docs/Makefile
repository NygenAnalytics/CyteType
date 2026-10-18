UVRUN ?= uv run --group docs

.PHONY: help html preview check-reference execute-docs

help:
	@echo 'make html                         Render saved notebook outputs and the API'
	@echo 'make preview                      Allow clearly labelled unexecuted pages'
	@echo 'make check-reference              Build and check the Sphinx API'
	@echo 'make execute-docs PAGES="quickstart"  Execute notebooks locally (real CyteType jobs)'

check-reference:
	@rm -rf build/api build/api-text
	@$(UVRUN) sphinx-build -b html -n -W source/reference/api build/api
	@$(UVRUN) python check_reference_coverage.py
	@$(UVRUN) sphinx-build -b text -n -W -d build/api/.doctrees source/reference/api build/api-text

html: check-reference
	@cd .. && $(UVRUN) python -m docs.build_docs

preview: check-reference
	@cd .. && $(UVRUN) python -m docs.build_docs --preview

execute-docs:
	$(if $(strip $(PAGES)),,$(error PAGES is required; name the notebooks to execute))
	@cd .. && $(UVRUN) python -m docs.execute_docs --pages "$(strip $(PAGES))"
