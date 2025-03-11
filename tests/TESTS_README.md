## Testing Differential Expression

The tests are in the `tests` directory and are run using Python's `pytest` library. This requires a one-time setup of the environment:

```bash
bash tests/setup_tests.sh
```

### Seting up the `renv` environment to run the tests

```bash
cd scripts/expression

module load rocker/rver/4.4.0
R 
# Then Run 
renv::init()
renv::restore()

```

After this, you can run the tests using from the main project directory: 

```bash
bash tests/run_tests.sh
```