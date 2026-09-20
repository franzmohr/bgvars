# bgvars

R package for the specification and estimation of global vector autoregressive models.

## Branching workflow

- `main` **is the development version**. It carries a development version number
  between releases, such as `1.2.3.9000`.
- New work goes on a short-lived branch created from `main` and is merged back
  into `main` **through a pull request**. Nothing is pushed to `main` directly,
  including small changes: the repository requires a pull request, and a direct
  push only succeeds by bypassing that rule.
- A release is a **tag plus a GitHub release**, not a branch and not a state of
  `main`. The released version is whatever the tag points at, so nothing is held
  back from `main` on the grounds that it is unreleased. A release candidate is
  the same thing, tagged `vX.Y.Z-rcN` and marked as a pre-release on GitHub.

### Cutting a release

On a branch off `main`, in one pull request:

1. `DESCRIPTION`: set `Version` to the release version, without the
   development suffix.
2. `cran-comments.md`: describe the submission that is actually being made,
   and confirm the check results against a fresh run rather than an older one.

Merge that pull request, then tag the merge commit and create the GitHub release
from the tag. The `release-version` workflow checks the tag against `DESCRIPTION`,
ignoring any `-rcN` suffix, so a mismatch fails the tag.

Immediately afterwards, open a second pull request setting `DESCRIPTION` back to
a development version, `X.Y.Z.9000`.
