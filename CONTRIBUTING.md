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


### Creating new releases

Creating a release and publishing a new heatrapy version into pypi is
only possible from the `master` branch. To do so, one needs to create a
new tag, relying on [Semantic Versioning](http://semver.org/), e.g.,
v2.0.9. There is a workflow that is triggered by pushing new tags into
the `master` branch. Merging a pull request into the `develop` branch
also triggers a publishing release, however in the test.pypi
repository. The publishing in test.pypi is only successful once per
tag version.


### Documentation
The [documentation](https://djsilva99.github.io/heatrapy) is generated
by the heatrapy docstrings and pushed in the branch `gh-pages`. You
can find how to generate the documentation and publish it in the
[README](https://github.com/djsilva99/heatrapy/tree/gh-pages) file of
that branch. There is also an additional branch for testing the
generated documentation: `gh-pages-dev`.
