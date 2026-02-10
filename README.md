# GEB Puzzles

This repo contains some code related to puzzles from Gödel, Escher, Bach, the nonfiction book by Doug Hofstadter. You can read a bit more about it [here](https://christian-johnson.github.io/posts/2026-02-09-geb/).

## Typogenetics
I set up a GitHub Actions workflow such that if a pull request is submitted that modifies `typogenetics/selfreps.txt`, an automated test will be run to determine whether the sequence submitted is self-replicating or not.

If you are interested in running the code itself to inspect a given sequence of interest, you use `uv` like this:
```
uv run typogenetics/main.py run-iterations --starting-string ATCGATCGATCG --verbose True --max-generations 5
```

Or you can loop through all possible sequences of a given length, you can do:
```
uv run typogenetics/main.py run-many-iterations --strand-length 10 --max-generations 3
```
