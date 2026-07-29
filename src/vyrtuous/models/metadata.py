def metadata(**kwargs):
    def decorator(func):
        for key, value in kwargs.items():
            setattr(func, key, value)
        return func

    return decorator


def command_metadata(**kwargs):
    def decorator(func):
        for key, value in kwargs.items():
            setattr(func, key, value)
        command = getattr(func, "__command__", None)
        if command:
            for key, value in kwargs.items():
                setattr(command, key, value)
        return func

    return decorator
