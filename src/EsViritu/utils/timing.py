import time
import functools
import logging

def timed_function(logger=None, level=logging.INFO):
    """
    Decorator to log the start, end, and execution time of a function.
    
    Args:
        logger: A logging.Logger instance. If None, uses the root logger.
        level: The logging level to use.
    
    Returns:
        The decorated function.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Get the appropriate logger
            log = logger or logging.getLogger("esv_logger")
            
            # Log start
            func_name = func.__name__
            log.log(level, f"Starting {func_name}...")
            
            # Track time
            start_time = time.time()
            
            # Call the function
            result = func(*args, **kwargs)
            
            # Calculate elapsed time
            elapsed_time = time.time() - start_time
            
            # Log completion with time
            log.log(level, f"Completed {func_name} in {elapsed_time:.2f} seconds")
            
            return result
        return wrapper
    return decorator
