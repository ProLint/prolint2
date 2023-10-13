import pytest
from bottle import HTTPResponse

from prolint2.server.server import ProLintDashboard


# Define a fixture for the app
@pytest.fixture
def app():
    return ProLintDashboard(debug_bool=True)


# Test the setup_routes() method
def test_setup_routes(app):
    # Call the method
    app.setup_routes()

    # Assert that the routes are correctly set up
    assert len(app.app.routes) == 20


# Test the serve_app() method
def test_serve_app(app):
    # Call the method
    response = app.serve_app()

    # Assert that the response is of type HTTPResponse
    assert isinstance(response, HTTPResponse)
    assert response.status_code == 200
    assert response.content_type == "text/html; charset=UTF-8"


# Add more test functions for other methods and functionalities
# TODO: Add more tests for the ProLintDashboard class
