""":mod:`cpg_islands.qt.composers` --- Functions to create Qt MVP triads
"""

from cpg_islands.models import (AppModel,
                                SeqInputModel,
                                ResultsModel,
                                EntrezModel)
from cpg_islands.qt.views import (AppView,
                                  SeqInputView,
                                  ResultsView,
                                  EntrezView)
from cpg_islands.presenters import (AppPresenter,
                                    SeqInputPresenter,
                                    ResultsPresenter,
                                    EntrezPresenter)


def create_app_presenter(argv):
    """Create a presenter with a Qt view.

    :return: the created presenter
    :rtype: :class:`AppPresenter`
    """
    results_model = ResultsModel()
    results_view = ResultsView()
    seq_input_model = SeqInputModel(results_model)
    seq_input_view = SeqInputView()
    entrez_model = EntrezModel(seq_input_model)
    entrez_view = EntrezView()
    app_model = AppModel(seq_input_model, entrez_model)
    app_view = AppView(entrez_view, seq_input_view, results_view)

    seq_input_presenter = SeqInputPresenter(seq_input_model, seq_input_view)
    results_presenter = ResultsPresenter(results_model, results_view)
    entrez_presenter = EntrezPresenter(entrez_model, entrez_view)
    app_presenter = AppPresenter(app_model, app_view)
    for object in [app_model, results_presenter,
                   seq_input_presenter, entrez_presenter,
                   app_presenter]:
        object.register_for_events()
    app_model.run(argv)

    return app_presenter
