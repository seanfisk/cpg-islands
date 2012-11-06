""":mod:`cpg_islands.qt.composers` --- Functions to create Qt MVP triads
"""

from cpg_islands.models import (ApplicationModel,
                                SequenceInputModel,
                                ResultsModel)
from cpg_islands.qt.views import (ApplicationView,
                                  SequenceInputView,
                                  ResultsView)
from cpg_islands.presenters import (ApplicationPresenter,
                                    SequenceInputPresenter,
                                    ResultsPresenter)


def create_application_presenter(argv):
    """Create a presenter with a Qt view.

    :return: the created presenter
    :rtype: :class:`ApplicationPresenter`
    """
    results_model = ResultsModel()
    results_view = ResultsView()
    results_presenter = ResultsPresenter(results_model, results_view)
    sequence_input_model = SequenceInputModel(results_model)
    sequence_input_view = SequenceInputView()
    sequence_input_presenter = SequenceInputPresenter(sequence_input_model,
                                                      sequence_input_view)
    application_model = ApplicationModel(sequence_input_model)
    application_view = ApplicationView(sequence_input_view, results_view)

    application_presenter = ApplicationPresenter(application_model,
                                                 application_view)
    application_model.run(argv)
    return application_presenter
