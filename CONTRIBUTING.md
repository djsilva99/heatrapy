# Contributing

Feel free to actively contribute to this project by opening pull
requests. They will all be considered. You can find suggested
contributions in the
[roadmap](https://github.com/djsilva99/heatrapy/wiki) by solving
[issues](https://github.com/djsilva99/heatrapy/issues).


### Development workflow

There are two main branches: `master`, which is the base branch, and
`develop`, which is the branch where new features are merged. All
created pull requests should be merged into the `develop` branch. From
time to time, the content of the `develop` branch is merged into the
`master` branch, where new releases are created.

To build the heatrapy package locally, first create a virtual
environment and install the dependencies:
```bash
$ python -m venv env
$ source env/bin/activate
(env) $ pip install -r requirements
```

Tests rely on pytest and can be executed with
```bash
(env) $ pytest
```

Currently, the coverage must exceed 70%, but this value will be
increased to 85% in a near future (according to the roadmap).

Lint relies on flake8, and can be validated locally with
```bash
(env) $ flake8 . --count --show-source --statistics
```

There is a workflow that validates the tests, coverage and lint before
merging pull requests. They must all pass.


### Documentation
The
[documentation](https://heatrapy.readthedocs.io/en/latest/index.html)
is generated automatically whenever new code is merged into the master
branch (which occurs every new release) by the read the docs
workflow. The content of the documentation takes into account the
information in the docs directory and the heatrapy docstrings. The
workflow was built using [this read the docs
tutorial](https://docs.readthedocs.io/en/stable/tutorial/index.html)


### Creating new releases

Creating a release and publishing a new heatrapy version into pypi is
only possible from the `master` branch. To do so, one needs to create
a new tag, relying on [Semantic Versioning](http://semver.org/), e.g.,
v2.0.9. There is a workflow that is triggered by pushing new tags into
the `master` branch. Merging a pull request into the `develop` branch
also triggers a publishing release, however in the test.pypi
repository. The publishing in test.pypi is only successful once per
tag version. Do not forget to change the new version tag in the
documentation before merging code into the master branch and creating
the version tag, namely in the following locations:

- setup.py
- README.md
- heatrapy./\_\_init\_\_.py
- .readthedocs.yaml
- docs/conf.py
- docs/index.rst

Thus, the steps that create a new release are the following:

1. Create a PR that adds the new version tag into the locations
   pointed out above. Merge it into the develop branch.
2. Merge into the develop branch into the master branch.
3. Create a tag with the new version in the master branch.
4. Test the install of the new version using pip.
