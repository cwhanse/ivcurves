import argparse
from pathlib import Path


DOCS_SUBMISSIONS = Path(__file__).parent / 'submissions'
newline_indented = '\n   '


def path_without_extension(file):
    return file.parent / file.stem


def create_py_file_rst(file):
    file = file.resolve()
    submissions_idx = file.parts.index('submissions')
    contents = f'''
{file.stem}
{'=' * len(file.stem)}

.. automodule:: {'.'.join(path_without_extension(file).parts[submissions_idx:])}
   :members:

'''.strip() + '\n\n'
    with open(f'{path_without_extension(file)}.rst', 'w') as f:
        f.write(contents)


def create_index_rst(file):
    file = file.resolve()
    is_one_subdir_down = lambda p: len(p.parts) >= 3 and p.parts[-3] == parent.stem
    while file.parts[-2] != 'submissions':
        parent = file.parent
        sub_idx_rst = (i for i in sorted(parent.rglob('index.rst')) if is_one_subdir_down(i))
        contents = f'''
{parent.stem}
{'=' * len(parent.stem)}

.. toctree::
   :maxdepth: 2

   {newline_indented.join(f.stem for f in sorted(parent.glob('*.rst')) if f.stem != 'index')}

   {newline_indented.join(f'{i.parent.stem}/index' for i in sub_idx_rst)}

'''.strip() + '\n\n'
        with open(parent / 'index.rst', 'w') as f:
            f.write(contents)
        file = parent


def add_to_submissions_index(gh_username):
    index_entry = f'{gh_username}/index'
    with open(DOCS_SUBMISSIONS / 'index.rst') as f:
        contents = f.read()
    if index_entry not in contents:
        with open(DOCS_SUBMISSIONS / 'index.rst', 'a') as f:
            f.write(f'   {index_entry}\n\n')


def tree(dir_path: Path, prefix: str=''):
    """A recursive generator, given a directory Path object
    will yield a visual tree structure line by line
    with each line prefixed by the same characters
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
                    'unless it is index.rst'
    )
    parser.add_argument('submission_path', type=Path,
                        help='The path to the directory containing the submission.')
    return parser


if __name__ == '__main__':
    args = get_argparser().parse_args()

    sp = args.submission_path

    if not sp.exists():
        raise ValueError('The given path does not exist.')
    if not sp.is_dir():
        raise ValueError('The given path must be a directory, not a file.')

    submissions_idx = sp.parts.index('submissions')
    gh_username = sp.parts[submissions_idx + 1]

    py_files = sorted(sp.rglob('*.py'))
    rst_files = []
    for f in py_files:
        parent = DOCS_SUBMISSIONS.joinpath(*f.parts[submissions_idx + 1:]).parent
        rst_files.append(parent / f'{f.stem}.rst')

    for f in rst_files:
        f.parent.mkdir(parents=True, exist_ok=True)
        # if not file.exists():
        create_py_file_rst(f)

    for f in rst_files:
        create_index_rst(f)

    add_to_submissions_index(gh_username)

    if len(rst_files) != 0:
        print('\n'.join(tree(DOCS_SUBMISSIONS / gh_username)))

