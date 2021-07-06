# Pull Requests

All pull requests should be reviewed by a core developer, and may include a review by a subject matter expert if the area of the PR is outside that of one of the core developers. In that case, the core developers will primarily review style and design, rather than substance.

Every PR to PowerModelsDistribution should strive to meet the following guidelines.

## PR Title

- Should be concise and clear, describing in a phrase the content of the PR
- Should include a prefix that describes the primary type of the PR
  - ADD: feature addition
  - FIX: bugfix
  - REF: refactor
  - UPD: updates to code for _e.g._ version bumps of dependencies
  - STY: style changes, no changes to function names, added features, etc.
  - DOC: documentation-only additions/changes
  - RM: dead code removal

## PR Body

- If the change is breaking, it should be clearly stated up front
- The purpose of this PR should be clearly stated right away
- Major changes / additions to the code should be summarized. In the case where a refactor was performed, the name changes of public functions should be documented in the body of the PR
- Any associated Issues should be referenced in the body of the PR, and it is accepted/encouraged to use Closes #XX to automatically close Issues after the PR is merged

## PR Code

- An entry should be added to CHANGELOG.md for every PR
- Documentation should be updated (See [Documentation](##Documentation) section above for guidelines)
- Unit tests should be added. In the case where existing unit tests were altered, an explanation for the change must be included
- Code should be rebased to the latest version of whatever branch the PR is aimed at (no merge conflicts!)

# Versions

PowerModelsDistribution follows the Semantic Versioning ([SemVer](https://semver.org/)) convention of `Major.minor.patch`, where `Major` indicates breaking changes, `minor` indicates non-breaking feature additions, and `patch` indicates non-breaking bugfixes.

Currently, because `Major==0`, `minor` indicates breaking changes and `patch` indicates any non-breaking change, including both feature additions and bugfixes. Once PowerModelsDistribution reaches `v1.0.0`, we will adhere strictly to the SemVer convention.

# Branch Management

The `main` branch is a [protected](https://help.github.com/en/github/administering-a-repository/about-protected-branches) branch, meaning that its history will always be contiguous and can never be overwritten.

Release candidate branches of the format `vM.m.0-rc` are also protected branches. These branches will contain only breaking changes and will not be merged into main until a new version is ready to be tagged. Pull requests including breaking changes should be directed into the next release candidate branch available, _e.g._ if the current version of the package is `v0.9.0`, the next release candidate branch will be `v0.10.0-rc`.

Pull requests that include only non-breaking changes can be merged directly into `main` once approved, and in the case of merge conflicts arising for release candidate branches, the `-rc` branch will need to be updated to include the latest `main`.

Pull requests will generally be merged using [squash and merge](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-request-merges#squash-and-merge-your-pull-request-commits) into the branch they are aimed at, with the exception of release candidate branches, which generally be merged using [rebase and merge](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-request-merges#rebase-and-merge-your-pull-request-commits) into main.
