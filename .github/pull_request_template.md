# Pull Request (PR) Template

Pull requests should be single-issue to the extent possible; by focusing a PR on a single issue/topic, code reviews are more manageable, new features are easier to track, code quality is easier to maintain, and technical debt is easier manage.

Every PR to PowerModelsDistribution should strive to meet the following guidelines.

## Title

- Should be concise and clear, describing in a phrase the content of the PR
- Should include a prefix that describes the primary type of the PR
  - ADD: feature addition
  - FIX: bugfix
  - REF: refactor
  - UPD: updates to code for _e.g._ version bumps of dependencies
  - STY: style changes, no changes to function names, added features, etc.
  - DOC: documentation-only additions/changes
  - RM: dead code removal

## Body

- If the change is breaking, it should be clearly stated up front
- The purpose of this PR should be clearly stated right away
- Major changes / additions to the code should be summarized. In the case where a refactor was performed, the name changes of public functions should be documented in the body of the PR
- Any associated Issues should be referenced in the body of the PR, and it is accepted/encouraged to use Closes #XX to automatically close Issues after the PR is merged

## Code

- Code changes should be relevant to the issue being addressed by the PR to the extent possible. Unrelated code changes are discouraged, and should be submitted in a separate PR
- An entry should be added to CHANGELOG.md for every PR
- Documentation should be updated (See CONTRIBUTING.md for guidelines)
- Unit tests should be added. In the case where existing unit tests were altered, an explanation for the change must be included
