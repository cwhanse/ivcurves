import argparse
from pathlib import Path


DOCS_SUBMISSIONS = Path(__file__).parent / 'submissions'


def path_without_extension(file):
    """
    Gets the entire file path, excluding the file extension.
    This is the complement of ``pathlib.Path.suffix``.

    Parameters
    ----------
    file: pathlib.Path
        The path.

    Returns
    -------
        pathlib.Path
    """
    return file.parent / file.stem


def create_py_file_rst(file):
    """
    Creates a .rst file to document a Python module using automodule.

    Parameters
    ----------
    file: pathlib.Path
        A path to the .rst file to be created.
    """
    submissions_idx = len(DOCS_SUBMISSIONS.parts) - 1
    contents = f"""
{file.stem}
{'=' * len(file.stem)}

.. automodule:: {'.'.join(path_without_extension(file).parts[submissions_idx:])}
   :members:

""".strip() + '\n\n'
    with open(f'{path_without_extension(file)}.rst', 'w') as f:
        f.write(contents)


def create_index_rst(file):
    """
    Creates an ``index.rst`` file in every subdirectory from
    ``DOCS_SUBMISSIONS`` to the .rst file pointed to by ``file``.

    Parameters
    ----------
    file: pathlib.Path
        A path pointing to a .rst file.
    """
    is_one_subdir_down = lambda p: len(p.parts) >= 3 and p.parts[-3] == parent.stem
    # backslashing are not allowed in f-strings
    newline_indented = '\n   '
    while file.parts[-2] != 'submissions':
        parent = file.parent
        sub_idx_rst = (i for i in sorted(parent.rglob('index.rst')) if is_one_subdir_down(i))
        contents = f"""
{parent.stem}
{'=' * len(parent.stem)}

.. toctree::
   :maxdepth: 2

   {newline_indented.join(f.stem for f in sorted(parent.glob('*.rst')) if f.stem != 'index')}

   {newline_indented.join(f'{i.parent.stem}/index' for i in sub_idx_rst)}

""".strip() + '\n\n'
        with open(parent / 'index.rst', 'w') as f:
            f.write(contents)
        file = parent


def add_to_submissions_index(gh_username):
    """
    Updates the ``submissions/index.rst`` file by adding the entry
    ``<gh_username>/index``, unless it already exists.
    """
    index_entry = f'{gh_username}/index'
    with open(DOCS_SUBMISSIONS / 'index.rst') as f:
        contents = f.read()
    if index_entry not in contents:
        with open(DOCS_SUBMISSIONS / 'index.rst', 'a') as f:
            f.write(f'   {index_entry}\n\n')


def tree(dir_path: Path, prefix: str=''):
    """
    A recursive generator, given a directory Path object
    will yield a visual tree structure line by line
    with each line prefixed by the same characters

    Parameters
    ----------
    dir_path: pathlib.Path
        A path.

    prefix: str, default ''
        The character prefix for each branch of the file tree.

    Returns
    -------
        generator(str)
    """
    # prefix components:
    space =  '    '
    branch = '│   '
    # pointers:
    tee =    '├── '
    last =   '└── '
    if prefix == '':
        yield dir_path.stem
    contents = sorted(dir_path.iterdir())
    # contents each get pointers that are ├── with a final └── :
    pointers = [tee] * (len(contents) - 1) + [last]
    for pointer, path in zip(pointers, contents):
        yield prefix + pointer + path.name
        if path.is_dir(): # extend the prefix and recurse:
            extension = branch if pointer == tee else space
            # i.e. space because last, └── , above so no more |
            yield from tree(path, prefix=prefix+extension)


def get_argparser():
    parser = argparse.ArgumentParser(
        description='Generates .rst files for documenting a submission.'
                    'Note a .rst file is not created if it already exists,'
                    'unless it is index.rst.'
    )
    parser.add_argument('pr_config_path', type=Path,
                        help='The path to the pr_config.json of a submission.')
    return parser


if __name__ == '__main__':
    args = get_argparser().parse_args()

    # resolve path to ensure 'submissions/<GitHub_username>' are in the path
    pr_cfg = args.pr_config_path.resolve()

    if not pr_cfg.exists():
        raise ValueError('The given path does not exist.')
    if pr_cfg.name != 'pr_config.json':
        raise ValueError('The given path must point to pr_config.json.')

    # submissions_gh = submissions/<GitHub_username>
    submissions_gh = pr_cfg.parent
    gh_username = submissions_gh.stem
    gh_username_idx = len(submissions_gh.parts) - 1

    py_files = sorted(submissions_gh.rglob('*.py'))
    # create a list of rst files to be created, excluding index.rst files
    rst_files = []
    for f in py_files:
        # parent is the path DOCS_SUBMISSIONS concatenated with the path to f
        # relative to the path submissions_gh
        parent = DOCS_SUBMISSIONS.joinpath(*f.parts[gh_username_idx:]).parent
        rst_files.append(parent / f'{f.stem}.rst')

    # create the rst files, if they do not already exist
    for f in rst_files:
        f.parent.mkdir(parents=True, exist_ok=True)
        if not f.exists():
            create_py_file_rst(f)

    # create the index.rst files
    for f in rst_files:
        create_index_rst(f)

    # add entry to 'submissions/index.rst'
    add_to_submissions_index(gh_username)

    if len(rst_files) != 0:
        print('\n'.join(tree(DOCS_SUBMISSIONS / gh_username)))

