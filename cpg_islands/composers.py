""":mod:`cpg_islands.composers` --- Functions to create MVP triads
"""

from models import ApplicationModel
from views.cli import CliApplicationView
from presenters import ApplicationPresenter


def create_cli_presenter():
    """Create a presenter with a command-line view.

    :return: the created presenter
    :rtype: :class:`ApplicationPresenter`
    """
    model = ApplicationModel()
    view = CliApplicationView()
    presenter = ApplicationPresenter(model, view)
    presenter.register_for_events()
    model.run()
    return presenter
