Changelog
=========

1.0.0 (2018-08-09)
------------------
- Initial release.

1.1.1 (2018-08-18)
------------------
- Bug-fix: When choosing the type of redundancy in the RedundancyType enum, TARGET_FLEXIBILITY and SOURCE_FELXIBILITY used to always equal to FLEXIBILITY, due to the scriptiness of Python. Now they are treated differently.
- Added partial redundancy to the the RedundancyType enum, plus logic for the Redundancy and RedundancyContribution classes.
- Added export of dictionaries containing Neofunctionalisation/FunctionChange/Element objects into a formatted HTML file.
- Added three more experiments involving the above.
