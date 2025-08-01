#!/usr/bin/env python3
"""
ProLint2 Test Runner
Comprehensive test execution with reporting and analysis.
"""

import sys
import os
import subprocess
import argparse
import time
from pathlib import Path
from typing import List, Dict, Optional

# Add prolint2 to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))


def run_command(cmd: List[str], capture_output: bool = True) -> tuple:
    """Run a command and return result."""
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd, 
            capture_output=capture_output, 
            text=True, 
            check=False
        )
        return result.returncode, result.stdout, result.stderr
    except Exception as e:
        print(f"Error running command: {e}")
        return 1, "", str(e)


def run_unit_tests(args) -> bool:
    """Run unit tests with pytest."""
    print("\n" + "="*50)
    print("RUNNING UNIT TESTS")
    print("="*50)
    
    cmd = [
        sys.executable, "-m", "pytest", 
        "prolint2/tests/",
        "-v",
        "--tb=short"
    ]
    
    if args.coverage:
        cmd.extend(["--cov=prolint2", "--cov-report=term-missing", "--cov-report=html"])
    
    if args.parallel:
        cmd.extend(["-n", "auto"])
    
    if args.markers:
        cmd.extend(["-m", args.markers])
    
    returncode, stdout, stderr = run_command(cmd, capture_output=False)
    
    if returncode == 0:
        print("✅ Unit tests PASSED")
        return True
    else:
        print("❌ Unit tests FAILED")
        return False


def run_integration_tests(args) -> bool:
    """Run integration tests."""
    print("\n" + "="*50)
    print("RUNNING INTEGRATION TESTS")
    print("="*50)
    
    cmd = [
        sys.executable, "-m", "pytest", 
        "prolint2/tests/",
        "-v",
        "-m", "integration",
        "--tb=short"
    ]
    
    returncode, stdout, stderr = run_command(cmd, capture_output=False)
    
    if returncode == 0:
        print("✅ Integration tests PASSED")
        return True
    else:
        print("❌ Integration tests FAILED")
        return False


def run_linting(args) -> bool:
    """Run code linting."""
    print("\n" + "="*50)
    print("RUNNING CODE LINTING")
    print("="*50)
    
    success = True
    
    # Flake8
    print("\n--- Running flake8 ---")
    cmd = ["flake8", "prolint2", "--max-line-length=88", "--extend-ignore=E203,W503"]
    returncode, stdout, stderr = run_command(cmd)
    if returncode != 0:
        print("❌ flake8 failed:")
        print(stdout)
        print(stderr)
        success = False
    else:
        print("✅ flake8 passed")
    
    # Black formatting check
    print("\n--- Checking black formatting ---")
    cmd = ["black", "--check", "--diff", "prolint2"]
    returncode, stdout, stderr = run_command(cmd)
    if returncode != 0:
        print("❌ Black formatting check failed:")
        print(stdout)
        success = False
    else:
        print("✅ Black formatting check passed")
    
    # isort import sorting check
    print("\n--- Checking isort import sorting ---")
    cmd = ["isort", "--check-only", "--diff", "prolint2"]
    returncode, stdout, stderr = run_command(cmd)
    if returncode != 0:
        print("❌ isort check failed:")
        print(stdout)
        success = False
    else:
        print("✅ isort check passed")
    
    return success


def run_type_checking(args) -> bool:
    """Run type checking with mypy."""
    print("\n" + "="*50)
    print("RUNNING TYPE CHECKING")
    print("="*50)
    
    cmd = ["mypy", "prolint2", "--ignore-missing-imports"]
    returncode, stdout, stderr = run_command(cmd)
    
    if returncode == 0:
        print("✅ Type checking PASSED")
        return True
    else:
        print("❌ Type checking FAILED:")
        print(stdout)
        print(stderr)
        return False


def run_security_check(args) -> bool:
    """Run security vulnerability check."""
    print("\n" + "="*50)
    print("RUNNING SECURITY CHECK")
    print("="*50)
    
    cmd = ["safety", "check"]
    returncode, stdout, stderr = run_command(cmd)
    
    if returncode == 0:
        print("✅ Security check PASSED")
        return True
    else:
        print("⚠️  Security check found issues:")
        print(stdout)
        print(stderr)
        return False  # Don't fail on security issues initially


def run_performance_tests(args) -> bool:
    """Run performance benchmarks."""
    print("\n" + "="*50)
    print("RUNNING PERFORMANCE TESTS")
    print("="*50)
    
    cmd = [
        sys.executable, "-m", "pytest", 
        "prolint2/tests/",
        "--benchmark-only",
        "-v"
    ]
    
    returncode, stdout, stderr = run_command(cmd, capture_output=False)
    
    if returncode == 0:
        print("✅ Performance tests PASSED")
        return True
    else:
        print("❌ Performance tests FAILED")
        return False


def run_examples(args) -> bool:
    """Test example scripts."""
    print("\n" + "="*50)
    print("TESTING EXAMPLE SCRIPTS")
    print("="*50)
    
    examples_dir = Path("examples")
    if not examples_dir.exists():
        print("⚠️  Examples directory not found")
        return True
    
    success = True
    example_files = [
        "logging_girk_demo.py",
        "validation_usage.py", 
        "memory_usage.py"
    ]
    
    for example in example_files:
        example_path = examples_dir / example
        if example_path.exists():
            print(f"\n--- Testing {example} ---")
            cmd = [sys.executable, str(example_path)]
            returncode, stdout, stderr = run_command(cmd)
            
            if returncode == 0:
                print(f"✅ {example} executed successfully")
            else:
                print(f"❌ {example} failed:")
                print(stderr)
                success = False
        else:
            print(f"⚠️  {example} not found")
    
    return success


def generate_report(results: Dict[str, bool], start_time: float):
    """Generate test report."""
    print("\n" + "="*60)
    print("TEST EXECUTION REPORT")
    print("="*60)
    
    total_tests = len(results)
    passed_tests = sum(1 for result in results.values() if result)
    failed_tests = total_tests - passed_tests
    
    print(f"Total test suites: {total_tests}")
    print(f"Passed: {passed_tests}")
    print(f"Failed: {failed_tests}")
    print(f"Success rate: {(passed_tests/total_tests)*100:.1f}%")
    print(f"Execution time: {time.time() - start_time:.2f} seconds")
    
    print("\nDetailed Results:")
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"  {test_name:<20} {status}")
    
    if failed_tests > 0:
        print(f"\n⚠️  {failed_tests} test suite(s) failed!")
        return False
    else:
        print("\n🎉 All test suites passed!")
        return True


def main():
    """Main test runner."""
    parser = argparse.ArgumentParser(description="ProLint2 Test Runner")
    
    parser.add_argument("--unit", action="store_true", help="Run unit tests")
    parser.add_argument("--integration", action="store_true", help="Run integration tests")
    parser.add_argument("--lint", action="store_true", help="Run linting")
    parser.add_argument("--type-check", action="store_true", help="Run type checking")
    parser.add_argument("--security", action="store_true", help="Run security check")
    parser.add_argument("--performance", action="store_true", help="Run performance tests")
    parser.add_argument("--examples", action="store_true", help="Test example scripts")
    parser.add_argument("--all", action="store_true", help="Run all tests")
    
    # Test configuration options
    parser.add_argument("--coverage", action="store_true", help="Generate coverage report")
    parser.add_argument("--parallel", action="store_true", help="Run tests in parallel")
    parser.add_argument("--markers", help="Pytest markers to select tests")
    
    # Output options
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--quick", action="store_true", help="Quick test run (essential only)")
    
    args = parser.parse_args()
    
    # If no specific tests selected, run essential tests
    if not any([args.unit, args.integration, args.lint, args.type_check, 
                args.security, args.performance, args.examples, args.all]):
        args.unit = True
        args.lint = True
    
    # If --all is specified, enable all tests
    if args.all:
        args.unit = True
        args.integration = True
        args.lint = True
        args.type_check = True
        args.security = True
        args.performance = True
        args.examples = True
    
    # If --quick is specified, run only essential tests
    if args.quick:
        args.unit = True
        args.lint = False
        args.type_check = False
        args.security = False
        args.performance = False
        args.examples = False
    
    start_time = time.time()
    results = {}
    
    print("ProLint2 Test Execution Started")
    print(f"Python: {sys.version}")
    print(f"Working directory: {os.getcwd()}")
    
    # Run selected test suites
    if args.unit:
        results["Unit Tests"] = run_unit_tests(args)
    
    if args.integration:
        results["Integration Tests"] = run_integration_tests(args)
    
    if args.lint:
        results["Linting"] = run_linting(args)
    
    if args.type_check:
        results["Type Checking"] = run_type_checking(args)
    
    if args.security:
        results["Security Check"] = run_security_check(args)
    
    if args.performance:
        results["Performance Tests"] = run_performance_tests(args)
    
    if args.examples:
        results["Example Scripts"] = run_examples(args)
    
    # Generate final report
    success = generate_report(results, start_time)
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
