from vyrtuous.utils.setup_logging import logger
import traceback
import time
import sys

# ANSI color codes for terminal output
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    DIM = '\033[2m'

class TestFailure(Exception):
    pass

class TestContext:
    def __init__(self, name):
        self.failures = []
        self.name = name
        self.passes = 0
        self.start = time.time()

    def assert_true(self, condition, msg):
        if not condition:
            raise TestFailure(msg)
        self.passes += 1

    def assert_equal(self, actual, expected, msg=None):
        if actual != expected:
            error_msg = msg or f"Expected {expected!r}, got {actual!r}"
            raise TestFailure(error_msg)
        self.passes += 1

    def fail(self, msg):
        raise TestFailure(msg)

    def record_exception(self, exc):
        tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        self.failures.append(tb)

    def summary(self):
        duration = time.time() - self.start
        return {
            "name": self.name,
            "passes": self.passes,
            "failures": len(self.failures),
            "duration": duration,
            "details": self.failures
        }

class TestReport:
    def __init__(self):
        self.errors = []
        self.failures = []
        self.passes = []
        self.start_time = time.time()
        self.current_section = None

    def start_section(self, section_name):
        """Start a new test section (like a test file in pytest)."""
        self.current_section = section_name
        print(f"\n{section_name} ", end='', flush=True)

    def error(self, name, exc):
        """Record an error (exception during test execution)."""
        self.errors.append((name, exc))
        print(f"{Colors.RED}E{Colors.RESET}", end='', flush=True)

    def fail(self, name, detail=None):
        """Record a failure (assertion failed)."""
        self.failures.append((name, detail))
        print(f"{Colors.RED}F{Colors.RESET}", end='', flush=True)

    def pass_(self, name):
        """Record a passing test."""
        self.passes.append(name)
        print(f"{Colors.GREEN}.{Colors.RESET}", end='', flush=True)

    def skip(self, name, reason=None):
        """Record a skipped test."""
        print(f"{Colors.YELLOW}s{Colors.RESET}", end='', flush=True)

    def print_summary(self):
        """Print pytest-style summary at the end."""
        print("\n")
        duration = time.time() - self.start_time
        
        # Print failures section
        if self.failures:
            print(f"\n{Colors.RED}{Colors.BOLD}{'='*70}")
            print("FAILURES")
            print(f"{'='*70}{Colors.RESET}")
            for i, (name, detail) in enumerate(self.failures, 1):
                print(f"\n{Colors.RED}{Colors.BOLD}_ {name} _{Colors.RESET}")
                if detail:
                    print(detail)
        
        # Print errors section
        if self.errors:
            print(f"\n{Colors.RED}{Colors.BOLD}{'='*70}")
            print("ERRORS")
            print(f"{'='*70}{Colors.RESET}")
            for i, (name, exc) in enumerate(self.errors, 1):
                print(f"\n{Colors.RED}{Colors.BOLD}_ {name} _{Colors.RESET}")
                tb = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
                print(tb)
        
        # Print summary line
        print(f"\n{Colors.BOLD}{'='*70}{Colors.RESET}")
        
        total = len(self.passes) + len(self.failures) + len(self.errors)
        parts = []
        
        if self.failures:
            parts.append(f"{Colors.RED}{len(self.failures)} failed{Colors.RESET}")
        if self.errors:
            parts.append(f"{Colors.RED}{len(self.errors)} error{Colors.RESET}")
        if self.passes:
            parts.append(f"{Colors.GREEN}{len(self.passes)} passed{Colors.RESET}")
        
        summary = ", ".join(parts)
        print(f"{summary} in {duration:.2f}s")
        
        # Return exit code (0 = success, 1 = failures)
        return 0 if not (self.failures or self.errors) else 1

    def print_short_summary(self):
        """Print a short summary (like pytest -v)."""
        print("\n")
        
        if self.passes:
            print(f"{Colors.GREEN}{Colors.BOLD}PASSED{Colors.RESET}")
            for name in self.passes:
                print(f"  {Colors.GREEN}✓{Colors.RESET} {name}")
        
        if self.failures:
            print(f"\n{Colors.RED}{Colors.BOLD}FAILED{Colors.RESET}")
            for name, detail in self.failures:
                print(f"  {Colors.RED}✗{Colors.RESET} {name}")
                if detail:
                    print(f"    {Colors.DIM}{detail}{Colors.RESET}")
        
        if self.errors:
            print(f"\n{Colors.RED}{Colors.BOLD}ERRORS{Colors.RESET}")
            for name, exc in self.errors:
                print(f"  {Colors.RED}E{Colors.RESET} {name}")
                print(f"    {Colors.DIM}{exc}{Colors.RESET}")
