from subprocess import call

def check_install():
    call("pytest $(ls sarabande/tests/*_tests.py)", shell=True)