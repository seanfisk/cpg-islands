class ApplicationPresenter(object):
    def __init__(self, model, view):
        self.model = model
        self.view = view

    def register_for_events(self):
        pass
