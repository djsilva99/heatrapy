# Change Log
All changes to this project are documented in this file, starting with version 2.0.0


## [2.0.3] - 2022-05-05

Fixed 2D latent heat computation for the 'explicit_k(x)' solver.

### Fixed
- 2D latent heat bug for the 'explicit_k(x)' solver.

### Changed
- default solver of SingleObject2D compute method is now 'implicit_k(x)'.
- updated vaccum thermal properties.


## [2.0.2] - 2022-04-13

Fixed latent heat and path bugs.

### Fixed
- Latent heat bug for 2D models
- Paths for windows
- Dependabot alerts


## [2.0.1] - 2021-03-29

Updated docstrings.

### Changed
- all docstrings to meet v2.0.0


## [2.0.0] - 2021-03-26

First version with two-dimensional features.

### Added
- SingleObject2D
- SystemObject2D

### Changed
- SingleObject1D
- SystemObject1D

### Fixed
- Data verification for all the classes
- Added docstrings in all methods
