# Copyright 2025, Drift+Noise GmbH

# This file is part of PRIIMA.
# PRIIMA is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# PRIIMA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with
# PRIIMA. If not, see https://www.gnu.org/licenses/gpl-3.0.html.

.PHONY: clean tags fast_cast fast_cast_debug lint

IS_IN_VIRTUALENV := $(shell python3 -c 'import sys; print(hasattr(sys, "real_prefix"))')
ifneq "$(IS_IN_VIRTUALENV)" "True"
    ACTIVATE_VENV := . venv/bin/activate;
endif

all: fast_cast test # build fast cast library first to have all imports ready

test: unittest

check: test lint

# default: do not build in debug mode!
fast_cast:
	mkdir -p fast_cast/build
	$(ACTIVATE_VENV) cd fast_cast/build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j

fast_cast_debug:
	mkdir -p fast_cast/build
	$(ACTIVATE_VENV) cd fast_cast/build && cmake -DCMAKE_BUILD_TYPE=Debug .. && make -j

unittest:
	$(ACTIVATE_VENV) pytest --color yes tests/*py

cover:  cover-clean
	$(ACTIVATE_VENV) pytest --cov=priima --color yes --cov-report=html tests/

libyear:
	$(ACTIVATE_VENV) libyear --sort -r requirements.txt

pylint:
	$(ACTIVATE_VENV) pylint priima/ tests/

hadolint:
	$(ACTIVATE_VENV) hadolint docker/priima/Dockerfile

check_style:
	$(ACTIVATE_VENV) flake8 priima/ tests/

lint: hadolint check_style pylint libyear

docker-image:
	docker-compose build priima

tags:
	ctags $$(find ./ -name '*.py' | grep -v venv)

clean:
	cd fast_cast/build && make clean

cover-clean:
	rm -rf htmlcov
	rm -f .coverage
