import warnings

warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in divide",
    category=RuntimeWarning
)
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in scalar divide",
    category=RuntimeWarning
)
warnings.filterwarnings(
    "ignore",
    message="divide by zero encountered in divide",
    category=RuntimeWarning
)