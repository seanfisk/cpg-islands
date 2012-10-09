Model
    The data and logic that serves the needs of a presenter for some aspect of the business.
    
    Business logic.
    
    Holds no reference to presenter.

View
    A portion of the interface of the application - what the user sees and interacts with.
    
    As thin as possible, basically just passthrough of events and getters and setters of display information.
    
    Holds no reference to presenter.
    
    Now has ability to do manual-only testing.
    
    Should only accept primitive types as parameters to function. Accepting complex types would require processing, and therefore tests.
    
Presenter
    Behavior that corresponds directly to customer stories. Application wiring.
    
    Holds references to model and view, possibly through constructor dependency injection.
    
    Interpreter of events.
    
    No public API.
    
    Stateless.
    
    Should not use any GUI toolkit types or methods.
