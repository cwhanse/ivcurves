import argparse
import csv
import json

from ivcurves import utils


def load_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)


def save_json(json_dict, filename):
    with open(filename, 'w') as file:
        return json.dump(json_dict, file, indent=2)


def load_overall_scores(filename):
    """
    Buildings a dictionary from test set filenames (excluding file extensions)
    to strings represented scores by reading from a CSV file with these
    columns: test_set, and score.

    Returns an empty dictionary if a FileNotFoundError is raised.

    Parameters
    ----------
    filename : str
        The path to the CSV file. The file extension must be included.

    Returns
    -------
    dict
        A dictionary and an empty string.
        If a FileNotFoundError was raised, an empty dictionary is returned.
    """
    overall_scores = {}
    try:
        with open(filename, newline='') as file:
            reader = csv.DictReader(file)
            for row in reader:
                overall_scores[row['test_set']] = row['score']
        return overall_scores
    except FileNotFoundError:
        return {}


def validate_overall_scores(overall_scores):
    """
    Given a dictionary of overall scores (see :func:`load_overall_scores`),
    the following are validated:

    - There is a score for every test set.
    - Each score is parsable as a float.

    Otherwise, a ``ValueError`` is raised.

    Parameters
    ----------
    overall_scores : dict
        A dictionary from test set filenames (excluding file extensions) to
        strings representing scores.

    Returns
    -------
    bool, str
        Boolean for whether ``overall_scores`` is valid or not, and a string
        for an error message.
    """
    overall_scores_keys = set(overall_scores.keys())

    # Proceed with normal validation
    valid_test_set_filenames = utils.get_filenames_in_directory(utils.TEST_SETS_DIR)
    missing_test_set_filenames = valid_test_set_filenames - overall_scores_keys

    if missing_test_set_filenames:
        names_list = list(missing_test_set_filenames)
        names_list.sort() # sort to keep same order in error message
        return False, f'Missing scores from these test sets: {", ".join(names_list)}'

    for name, score_str in overall_scores.items():
        try:
            float(score_str) # validate is a number
        except ValueError:
            return False, f"The score of test set '{name}' must parse to a float: {score_str}"

    return True, ''


def write_overall_scores_to_database(database, pr_number, pr_author,
                                     pr_closed_datetime, merge_commit,
                                     submission_main, overall_scores):
    """
    Writes an entry in the JSON scores database.
    The scores database has this structure:

    .. code-block::

       {
          "<pr_number>": {
              "username": "<pr_author>",
              "submission_datetime": "<pr_closed_datetime>",
              "merge_commit": "<pr_merge_commit_sha>",
              "submission_main": "<SUBMISSION_MAIN_from_pr_config.json>",
              "test_sets": {
                  "<test_set_filename>": "<score>",
                  ...
              }
          }
       }

    As tables, the database has two tables: submissions, and test_set_scores.
    The submissions table has these columns: ``pr_number``, ``username``,
    ``submission_datetime``, and ``test_sets``.
    The primary key is ``pr_number``.
    The test_set_scores table has these columns: ``test_set``, and ``score``.
    The primary key is a compound key of ``pr_number`` and ``test_set``.

    Parameters
    ----------
    database : dict
        The scores database.

    pr_number : int
        The pull request number.

    pr_author : str
        The author of the pull request.

    pr_closed_datetime : str
        The datetime the pull request was closed.

    merge_commit : str
        The SHA of the commit when the pull request was merged.

    submission_main : str
        SUBMISSION_MAIN from the submission's pr_config.json.

    overall_scores : dict
        A dictionary from test set filenames (excluding file extensions) to
        strings representing scores.
    """
    # only want strings as dict keys
    database[str(pr_number)] = {
        'username': pr_author,
        'submission_datetime': pr_closed_datetime,
        'merge_commit': merge_commit,
        'submission_main': submission_main,
        'test_sets': overall_scores
    }


def mark_submission_broken(database, pr_number, validation_msg):
    """
    Adds the entry ``"broken": true`` to indicate that the submission does
    not run will all the test sets.

    See :func:`write_overall_scores_to_database` for more info on ``database``.

    Parameters
    ----------
    database : dict
        The scores database.

    pr_number : int
        The pull request number.
    """
    database[str(pr_number)]['broken'] = validation_msg


def get_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pr-author', dest='pr_author', type=str,
                        help='GitHub username of the pull request author.')
    parser.add_argument('--pr-number', dest='pr_number', type=int,
                        help='GitHub pull request number.')
    parser.add_argument('--pr-closed-at', dest='pr_closed_datetime', type=str,
                        help='Datetime when the GitHub pull request closed.')
    parser.add_argument('--merge-commit', dest='merge_commit',
                        help='The commit created when the pull request was merged.')
    parser.add_argument('--submission-main', dest='submission_main',
                        help="The SUBMISSION_MAIN of the submission's pr_config.json")
    parser.add_argument('--overall-scores-path', dest='overall_scores_path',
                        type=str, help='Path to the CSV of overall scores.')
    parser.add_argument('--database-path', dest='database_path', type=str,
                        help='Path to the JSON scores database.')
    parser.add_argument('--broken-if-invalid', action=argparse.BooleanOptionalAction,
                        default=False,
                        help='Mark submission broken instead of raising an exception.')
    parser.add_argument('--save-database', action=argparse.BooleanOptionalAction,
                        default=True, help='Save updates to the database.')
    return parser


def record_scores():
    args = get_argparser().parse_args()

    overall_scores = load_overall_scores(args.overall_scores_path)
    has_valid_scores, validation_msg = validate_overall_scores(overall_scores)

    if not args.save_database:
        return

    database = load_json(args.database_path)

    if has_valid_scores:
        write_overall_scores_to_database(
            database, args.pr_number, args.pr_author, args.pr_closed_datetime,
            args.merge_commit, args.submission_main, overall_scores
        )
    elif args.broken_if_invalid:
        mark_submission_broken(database, args.pr_number, validation_msg)
    else:
        raise ValueError(validation_msg)

    save_json(database, utils.REPO_ROOT_DIR / args.database_path)


if __name__ == '__main__':
    record_scores()
