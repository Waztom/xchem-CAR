"""
Custom exceptions for the Recipe Generator module.
"""


class RecipeGeneratorError(Exception):
    """Base exception for recipe generator errors."""
    pass


class TemplateValidationError(RecipeGeneratorError):
    """Raised when a recipe template is invalid."""
    pass


class DesignMatrixError(RecipeGeneratorError):
    """Raised when a design matrix is invalid or cannot be parsed."""
    pass


class ActionNotFoundError(RecipeGeneratorError):
    """Raised when a referenced action cannot be found in the base recipe."""
    
    def __init__(self, action_id: str, session: int = None, actionnumber: int = None):
        self.action_id = action_id
        self.session = session
        self.actionnumber = actionnumber
        
        if session and actionnumber:
            msg = f"Action '{action_id}' not found at session {session}, actionnumber {actionnumber}"
        else:
            msg = f"Action '{action_id}' not found in recipe"
        super().__init__(msg)


class SessionNotFoundError(RecipeGeneratorError):
    """Raised when a referenced session cannot be found in the base recipe."""
    
    def __init__(self, session_id: str, sessionnumber: int = None):
        self.session_id = session_id
        self.sessionnumber = sessionnumber
        
        if sessionnumber:
            msg = f"Session '{session_id}' not found at sessionnumber {sessionnumber}"
        else:
            msg = f"Session '{session_id}' not found in recipe"
        super().__init__(msg)


class PathNotFoundError(RecipeGeneratorError):
    """Raised when a path within an action cannot be resolved."""
    
    def __init__(self, path: str, action_id: str):
        self.path = path
        self.action_id = action_id
        super().__init__(f"Path '{path}' not found in action '{action_id}'")


class OrderingError(RecipeGeneratorError):
    """Raised when action or session ordering is invalid."""
    pass
