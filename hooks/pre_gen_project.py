#!/usr/bin/env python3
"""Validate cookiecutter answers before the project tree is rendered."""

import re
import sys

project_slug = "{{ cookiecutter.project_slug }}"
contact_email = "{{ cookiecutter.contact_email }}"
first_system_name = "{{ cookiecutter.first_system_name }}"

errors = []

if not re.match(r"^[A-Za-z][A-Za-z0-9_\-]*$", project_slug):
    errors.append(
        f"project_slug '{project_slug}' is not valid: use letters, digits, "
        "'_' or '-', starting with a letter."
    )

if not re.match(r"^[A-Za-z][A-Za-z0-9_]*$", first_system_name):
    errors.append(
        f"first_system_name '{first_system_name}' is not valid: this becomes "
        "both a configs/<name>.yaml filename and a Python-safe identifier "
        "used elsewhere — use letters, digits or underscores, starting with "
        "a letter."
    )

if contact_email.strip() in ("", "you@example.org") or "@" not in contact_email:
    errors.append(
        "contact_email is required and must be a real address: EBI's "
        "InterProScan REST service (stage 1 annotation) rejects requests "
        "without one."
    )

if errors:
    print("\nERROR: fix the following before generating this project:\n")
    for e in errors:
        print(f"  - {e}")
    print()
    sys.exit(1)
