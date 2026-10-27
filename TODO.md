# TODO

Running list for the Fall 2026 (`f26`) revision. Check items off as they
land; add new ones at the bottom of the relevant section.

## Course framing

- [ ] Title page (`\titleTH` in `format106.sty`): says "NE 860 /
      Numerical Neutron Transport"; update to NE 806 / Neutronics. Drop
      the unused `\plogo` MIT logo and revisit "Created Spring 2011".
- [ ] Rewrite `introduction.tex` for NE 806 at K-State (currently the
      MIT 22.106 origin story).
- [ ] Remove references to the old course: "22.106 Stellar site"
      (`topics/criticality_safety.tex`), "in your homework"
      (`monte_carlo/tallies_and_statistics.tex`), MIT 18.335
      (`deterministic/discrete_ordinates.tex`).

## Part I: Nuclear data

- [ ] Resolve the transcription markup chapter by chapter, then delete
      it: 96 `\sourcepage`, 6 `\editorial`, 1 `\uncertain`, 8
      `\marginnoteinline`.
  - [ ] `neutron_scattering_1`
  - [ ] `quantum_scattering`
  - [ ] `neutron_scattering_2`
  - [ ] `reaction_theory_1`
  - [ ] `reaction_theory_2`
  - [ ] `quantum_angular_momentum`
  - [ ] `optical_model`
  - [ ] `r_matrix_theory`
  - [ ] `resonances_and_multiple_channels`
- [ ] Add `\label{lec:<file name>}` to each chapter (only
      `quantum_scattering` and `scattering_kinematics` have one).
- [ ] Add Blatt & Weisskopf to the bibliography citations; the notes
      say "mostly following B&W" without a `\cite`.
- [ ] Consider chapters on S(α,β)/thermal scattering and
      thermalization (old stubs were removed).

## Part II: Deterministic methods

- [ ] `multigroup_method`: body is a single paragraph; "To be
      continued" in the body and in Further Reading.
- [ ] `cpm_cylinder`: Sample Code is "To be continued"; Further Reading
      names Lewis & Miller, Hébert, and Carlvik (1964) without citing
      them.
- [ ] `pn_method_and_diffusion`: "A Simple Code" and Further Reading
      are "To be continued".
- [ ] `neutron_slowing_down`: Further Reading is "To be continued".
- [ ] Consider chapters on the method of characteristics, linearity and
      reciprocity, and coupling to other physics (old stubs were
      removed; git history has their one-line intros).

## Part III: Monte Carlo

- [ ] Merge the statistics material in `monte_carlo_introduction`
      (population moments) and `tallies_and_statistics` (sample
      statistics, CLT, batches) into one coherent treatment.
- [ ] Convert slide-derived text to prose: `tallies_and_statistics` and
      `random_numbers_and_collision_physics` still use `mitred` bullets
      and bare `{Heading}` lines left over from Beamer frames.
- [ ] `monte_carlo_introduction`: the Metropolis and Ulam photo is
      commented out and `images/monte_carlo/metropolis_and_ulam.jpg`
      does not exist; find the image or drop the figure.
- [ ] `automated_variance_reduction`: reads as a paper ("this paper");
      rework as a lecture. Add the six missing bib entries
      (`booth1985mcv`, `larsen2007hmc`, `peplow2009sot`,
      `peplow2010hmc`, `vanriper1997ava`, `wagner2009fwc`) and fix the
      overfull results tables (~60pt).
- [ ] Consider a Monte Carlo criticality (k-eigenvalue) chapter; the old
      `criticality.tex` stub held one exercise worth keeping.

## Part IV: Additional topics

- [ ] `criticality_safety`: overfull box (~36pt).
- [ ] Constructive solid geometry was a stub; fold into
      `computational_geometry` if wanted.

## Appendices

- [ ] Nomenclature: verify the † entries (Bell & Glasstone) against
      the editions on hand, then remove the daggers and the editorial
      note.
- [ ] Nomenclature: keep current as notation settles during revision;
      reconcile the notes' mixed use of `\mathbf{\Omega}` and
      `\hat{\Omega}`.

## Tooling and housekeeping

- [ ] `format106.sty` loads the obsolete `SIunits` and `times`
      packages; move to `siunitx` and `newtx`/`mathptmx`, and drop the
      `\renewcommand{\unit}` workaround in `preamble.tex`.
- [ ] `algorithm2e.sty` is vendored because the local TeX Live lacks
      it; drop it once installs include `texlive-science`.
- [ ] The criticality safety deck moved to
      `../slides/criticality_safety/`; commit it in the course repo. It
      still needs `mit-beamer-template` to build.
- [ ] Optional: CI that runs both `latexmk` builds and publishes the
      PDFs.
