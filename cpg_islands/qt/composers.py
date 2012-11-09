""":mod:`cpg_islands.qt.composers` --- Functions to create Qt MVP triads
"""

from cpg_islands.models import (AppModel,
                                SeqInputModel,
                                ResultsModel)
from cpg_islands.qt.views import (AppView,
                                  SeqInputView,
                                  ResultsView)
from cpg_islands.presenters import (AppPresenter,
                                    SeqInputPresenter,
                                    ResultsPresenter)


def create_app_presenter(argv):
    """Create a presenter with a Qt view.

    :return: the created presenter
    :rtype: :class:`AppPresenter`
    """
    results_model = ResultsModel()
    results_view = ResultsView()
    seq_input_model = SeqInputModel(results_model)
    seq_input_view = SeqInputView()
    app_model = AppModel(seq_input_model, results_model)
    app_view = AppView(seq_input_view, results_view)

    seq_input_presenter = SeqInputPresenter(seq_input_model, seq_input_view)
    results_presenter = ResultsPresenter(results_model, results_view)
    app_presenter = AppPresenter(app_model, app_view)

    for presenter in [app_model, results_presenter,
                      seq_input_presenter, app_presenter]:
        presenter.register_for_events()
    app_model.run(argv)

    return app_presenter
