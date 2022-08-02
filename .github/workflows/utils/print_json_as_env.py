import argparse
import json
import pathlib


def load_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)


def validate_pr_config(pr_config_json):
    """
    Validates that ``RUN_SCORER`` is a Boolean, and ``REQUIREMENTS`` and
    ``SUBMISSION_MAIN`` are paths that point to existing files.

    This can throw a ``ValueError``, ``FileNotFoundError``,
    or ``RuntimeError``.

    Parameters
    ----------
    pr_config_json : dict
        A dictionary from environment variable names to their values.

    Returns
    -------
    dict
        A validated mapping from environment variables to to their values.
    """
    pr_config_validated = {}
    valid_keys_to_value_types = {'RUN_SCORER': bool,
                                 'REQUIREMENTS': pathlib.Path,
                                 'SUBMISSION_MAIN': pathlib.Path}

    for key, value_type in valid_keys_to_value_types.items():
        pr_config_validated[key] = value_type(pr_config_json[key])

    for k, v in valid_keys_to_value_types.items():
        if isinstance(v, pathlib.Path):
            valid_keys_to_value_types[k] = f"'{v}'"
            valid_keys_to_value_types[f'{k}_PATH'] = f"'{v.parent}'"

    return pr_config_validated


def format_bool_variables(key, validated_dict):
    value = validated_dict[key]
    if isinstance(value, bool):
        validated_dict[key] = 'true' if value else 'false'


def format_path_variables(key, validated_dict, split_path_variables):
    quote_path = lambda path: f"'{path}'"
    value = validated_dict[key]
    if isinstance(value, pathlib.Path):
        validated_dict[key] = quote_path(value)
        if split_path_variables:
            validated_dict[f'{key}_FILENAME'] = quote_path(value.name)
            validated_dict[f'{key}_PATH'] = quote_path(value.parent)


def format_variable_values(validated_dict, split_path_variables=False):
    for k in list(validated_dict.keys()):
        format_bool_variables(k, validated_dict)
        format_path_variables(k, validated_dict, split_path_variables)


def print_json_as_env(validated_dict):
    for k, v in validated_dict.items():
        print(f'{k}={v}')


def get_argparser():
    parser = argparse.ArgumentParser(
        description='Prints entries in flat JSON object like environment variables.'
    )
    parser.add_argument('path', type=str, help='Path to the JSON file')
    parser.add_argument('--validate-pr-config', action=argparse.BooleanOptionalAction,
                        help='Runs the pr_config.json validator.')
    parser.add_argument('--split-path-variables', action=argparse.BooleanOptionalAction,
                        help='Adds two additional variables when a path variable P is encountered: a parent directory variable (P_PATH), and a filename variable (P_FILENAME).')
    return parser


if __name__ == '__main__':
    args = get_argparser().parse_args()
    flat_json = load_json(args.path)

    if args.validate_pr_config:
        flat_json = validate_pr_config(flat_json)

    # Run bash environment variable formatting
    format_variable_values(flat_json, split_path_variables=args.split_path_variables)

    print_json_as_env(flat_json)

