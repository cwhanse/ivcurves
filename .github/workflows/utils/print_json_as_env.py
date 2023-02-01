import argparse
from pathlib import Path

import ivcurves.utils as utils


def validate_pr_config(pr_config_json):
    """
    Validates that ``RUN_SCORER`` is a Boolean, and ``REQUIREMENTS`` and
    ``SUBMISSION_MAIN`` are paths that point to existing files, and converts
    the values of ``REQUIREMENTS`` and ``SUBMISSION_MAIN`` to ``pathlib.Path``.
    A ``ValueError`` is raised if the validation fails.

    Parameters
    ----------
    pr_config_json : dict
        A dictionary from environment variable names to their values.

    Returns
    -------
    dict
        A validated mapping from environment variables to to their values.
    """
    defn = lambda value_type, validator, error_msg: {
        'value_type': value_type,
        'validator': validator,
        'error_msg': error_msg
    }
    schema = {
        'RUN_SCORER': defn(bool, lambda v: isinstance(v, bool), 'must be a Boolean'),
        'REQUIREMENTS': defn(Path, lambda v: Path(v).exists(), 'path does not exist'),
        'SUBMISSION_MAIN': defn(Path, lambda v: Path(v).exists(), 'path does not exist')
    }

    missing_keys = set(schema.keys()) - set(pr_config_json.keys())
    if missing_keys:
        raise ValueError(f'Missing required keys in pr_config: {missing_keys}')
    extra_keys = set(pr_config_json.keys()) - set(schema.keys())
    if extra_keys:
        raise ValueError(f'Unknown keys in pr_config: {extra_keys}')

    for key, defn in schema.items():
        value_type, validator, error_msg = (
            defn['value_type'], defn['validator'], defn['error_msg']
        )
        if not validator(pr_config_json[key]):
            raise ValueError(f'{key}: {error_msg}')

        pr_config_json[key] = value_type(pr_config_json[key])


def format_bool_variables(key, validated_dict, **options):
    """
    Converts ``validated_dict[key]`` from a bool to a bash-style Boolean string
    If ``validated_dict[key]`` is not a bool, nothing is done. This modifies
    ``validated_dict``.

    Parameters
    ----------
    key : str
        A key of ``validated_dict``.

    validated_dict : dict
        A dict with key ``key`` where ``validated_dict[key]`` is possibly
        a bool.

    options : kwargs
        Any keyword arguments not used by this function.
    """
    value = validated_dict[key]
    if isinstance(value, bool):
        validated_dict[key] = 'true' if value else 'false'


def format_path_variables(key, validated_dict, split_path_variables=True, **options):
    """
    Converts ``validated_dict[key]`` from a ``pathlib.Path`` to a string. If
    ``validated_dict[key]`` is not a ``pathlib.Path``, nothing is done. This
    modifies ``validated_dict``.

    Parameters
    ----------
    key : str
        A key of ``validated_dict``.

    validated_dict : dict
        A dict with key ``key`` where ``validated_dict[key]`` is possibly
        a ``pathlib.Path``.

    split_path_variables : bool
        Adds two additional entires to ``validated_dict``:

            - ``{key}_FILENAME``: only the filename that
                ``validated_dict[key]`` points to.
            - ``{key}_PATH``: only the path to the parent folder containing the
                file ``validated_dict[key]`` points to.

    options : kwargs
        Any keyword arguments not used by this function.
    """
    value = validated_dict[key]
    if isinstance(value, Path):
        validated_dict[key] = str(value)
        if split_path_variables:
            validated_dict[f'{key}_FILENAME'] = str(value.name)
            validated_dict[f'{key}_PATH'] = str(value.parent)


def format_variable_values(validated_dict, **options):
    """
    Iterates through the keys of ``validated_dict`` and runs functions to
    modify their corresponding values. It allows for new keys to be added to
    ``validated_dict`` by the functions, but any new keys will not be iterated
    over.

    The functions called must have these three positional parameters:

    #. ``key``: The current key of the iteration.
    #. ``validated_dict``
    #. ``options``

    Parameters
    ----------
    validated_dict : dict
        A dict.

    options : kwargs
        A dict containing options that affect the behavior of the functions
        called during the iteration.
    """
    functions = [
        format_bool_variables,
        format_path_variables
    ]
    for k in list(validated_dict.keys()):
        for f in functions:
            f(k, validated_dict, **options)


def print_json_as_env(validated_dict):
    """
    Iterates through ``validated_dict`` and prints the string
    ``'{key}={value}'``.
    """
    for k, v in validated_dict.items():
        print(f'{k}={v}')


def get_argparser():
    parser = argparse.ArgumentParser(
        description='Prints entries in flat JSON object like environment variables.'
    )
    parser.add_argument('path', type=str, help='Path to the JSON file.')
    parser.add_argument('--validate-pr-config', action=argparse.BooleanOptionalAction,
                        help='Runs the pr_config.json validator.')
    parser.add_argument('--split-path-variables', action=argparse.BooleanOptionalAction,
                        help='Adds two additional variables when a path variable P is encountered: a parent directory variable (P_PATH), and a filename variable (P_FILENAME).')
    parser.add_argument('--quote-path-variables', action=argparse.BooleanOptionalAction,
                        help='Wrap path variables in double quotes (").')
    return parser


if __name__ == '__main__':
    args = get_argparser().parse_args()
    flat_json = utils.load_json(args.path)

    if args.validate_pr_config:
        validate_pr_config(flat_json)

    # Run bash environment variable formatting
    format_variable_values(flat_json, **args.__dict__)

    print_json_as_env(flat_json)
