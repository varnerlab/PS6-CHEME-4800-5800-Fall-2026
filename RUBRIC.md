# PS6 grading rubric

PS6 is worth a maximum of **4 points**: Parts 1–2 and the two questions in
[responses.md](responses.md). There are **40 tests**, 20 each in
[testme_part_1.jl](testme_part_1.jl) and [testme_part_2.jl](testme_part_2.jl).

| Score | Requirements |
|:---:|---|
| 0 | Nothing runs; no tests pass. |
| 1 | Some tests pass, but not a majority (1–20). |
| 2 | A majority of tests pass, but at least one fails (21–39). |
| 3 | All 40 tests pass, but questions 1 and 2 are not both satisfactorily answered. |
| 4 | All 40 tests pass and questions 1 and 2 are satisfactorily answered. |

Run `check_submission.jl` to see your test results and find unanswered
questions. Passing all the tests and the checker does not by itself earn a 4; the
instructor reads and grades your answers. An answer is satisfactory when it
answers every part of the question, gives the numbers it uses, and its reasoning
is correct.

To grade, the instructor copies your `src/Compute.jl`, any helper files, and
`Include.jl` into a clean copy of the assignment, reruns the tests, and then
reads your `responses.md`.

The policies in [README.md](README.md) also apply: you score `0` if you make no
initial submission or if you copy the reference solution.
