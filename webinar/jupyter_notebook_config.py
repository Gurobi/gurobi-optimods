# Configuration file for jupyter-notebook.

c = get_config()  # noqa

# Don't open the browser, notebook will be accessed via embedding
c.NotebookApp.open_browser = False

# Fix the port to 8888, don't allow retries (8888 is hardcoded in index.html)
c.NotebookApp.port = 8888
c.NotebookApp.port_retries = 0

# The IP address the notebook server will listen on.
c.NotebookApp.ip = "localhost"

# Operating in the background, don't ask questions
c.JupyterApp.answer_yes = True

# Disable all authentication and cross-site forgery checks to make
# embedding from localhost easy
c.NotebookApp.token = ""
c.NotebookApp.password_required = False
c.NotebookApp.disable_check_xsrf = True
c.NotebookApp.tornado_settings = {
    "headers": {
        "Content-Security-Policy": "frame-ancestors http://localhost:8000 'self' "
    }
}
