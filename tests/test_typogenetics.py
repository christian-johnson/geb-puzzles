from typogenetics.main import run_a_set_of_generations


def test_selfreps():
    with open("typogenetics/selfreps.txt") as f:
        selfreps = f.readlines()
    # strip newlines
    selfreps = [s.strip("\n") for s in selfreps]

    for selfrep in selfreps:
        assert run_a_set_of_generations(
            starting_string=selfrep, max_generations=3, verbose=False
        )
