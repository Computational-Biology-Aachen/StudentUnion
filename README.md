# StudentUnion
A Repository of all internim code shared between the students of the Computational Life Science Group at the RWTH Aachen

## Code Reviewing

To ensure a way to create code that actually is readable by other people, we want to introduce a code reviewing process between the students. Each student will have to upload code each week, that they were working on and let another person review it. That person will then have to be able to understand what the code does and be able to use it. They may find errors and have other recommendations of bettering the code. To make this process fun, the reviewer will get a randomized Bingo Card and if they can get a Bingo, the Author of the code will have to bring cake for the entire group.

### How to upload your code

To make the review as simple as possible, we will use the Pull Request method already integrated into Git. This will give the students also a small insight on how to use git pulls, pushes, branches and overall the group working idea behind a GitHub repository.

1. You should already find a branch with your name on it. Please clone the repository of that branch and only work on that. Do not push to main. (take care that the branch is up to date to `main`)
2. There should also be a folder with your name on it. You are free to upload your code there as you wish. You can format your directory as you want, the only rule is that you add a README.md to explain your directory and what to find where.
3. Once you are happy with your changes to your directory/branch, you should commit a pull-request to the `main` branch.

Here are the main things should include and have.

1. A clear README on what the code does, how to run it, and what files are included. You can choose if you do a new README for each new code, or you do one singular one. It should be clear though.
2. If your code includes functions, the names should be clear and they should include [MyPy type hints](https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html) and a [docstring](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings). There are several "templates" for docstring. We do not oblige you to follow one, however it should be clear, what the function does, needs, returns, and raises (if applicable). 
3. Comments are useful for more complex code. This is helpful for the reviewer and anybody else that wants to use/look at your code. It is also helpful when trying to tie the real world and what your code calculates.

```python
# WRONG
def func(
    x = ["ACCTGC", "ATTG", "CCTTA"]
):
    y = []
    z = []

    for i in range(len(x)):
        if len(x[i]) >= 5:
            y.append(x[i])
        else:
            z.append(x[i])

    return y

```

```python
# CORRECT
def sort_dna_correct_length(
    lst_dna: list[str] = ["ACCTGC", "ATTG", "CCTTA"],
    min_length: int = 5
):
    """Extracts DNA strands that are of a certain minimum length.
    
    Uses a list of DNA strands and a minimum length to extract only those strands that are of the minimum length or longer.
    
    Args:
        lst_dna (list[str], optional): List of DNA strands. Defaults to ["ACCTGC", "ATTG", "CCTTA"].
        min_length (int, optional): Minimum length of DNA strands to extract. Defaults to 5.
        
    Returns:
        list[str]: List of DNA strands that are of the minimum length or longer.
    """

    dna_correct_length = [] # List of correct lengths of DNA
    dna_wrong_length = [] # List of wrong lengths of DNA

    for idx in range(len(lst_dna)): # Iterate through list of DNA strands and extract ones that are the minimum length and put them into a new list
        if len(lst_dna[idx]) >= min_length:
            dna_correct_length.append(lst_dna[idx])
        else:
            dna_wrong_length.append(lst_dna[idx])

    return dna_correct_length

# EVEN BETTER
def sort_dna_correct_length(
    lst_dna: list = ["ACCTGC", "ATTG", "CCTTA"],
    min_length: int = 5
):
    """Extracts DNA strands that are of a certain minimum length.
    
    Uses a list of DNA strands and a minimum length to extract only those strands that are of the minimum length or longer.
    
    Args:
        lst_dna (list[str], optional): List of DNA strands. Defaults to ["ACCTGC", "ATTG", "CCTTA"].
        min_length (int, optional): Minimum length of DNA strands to extract. Defaults to 5.
        
    Returns:
        list[str]: List of DNA strands that are of the minimum length or longer.
    """
    return [dna for dna in lst_dna if len(dna) >= min_length]

```

## How to review your code

Each week you will be assigned to one of the students branches and their pull-request. You job is now to review their code and how they structured their directory. Here are things you need to keep an eye out for.

1. You can run the code.
   1. Every package is included in the `requirements` file.
   2. The README has a detailed walkthrough on how to run the code
2. Does the code include a significant README?
   1. What does the project do?
   2. What files are included?
   3. is there a walkthrough on how to get the code working?
3. Each function uses clever type hinting
4. Variables are understood from the names themselves
5. Can you rewrite the code to be cleaner?
   
To make the reviewing process more fun, we have created a bingo card, that is randomized for each reviewer. If you get a bingo (5 in a row horizontally, vertically or diagonally) the author of the code will have to bring cake next week.

### Bingo

To get Bingo to run, you have to have [npm](https://docs.npmjs.com/downloading-and-installing-node-js-and-npm) installed. You then have to run the following code inside of this directory:

```bash
npm install
npm run dev
```

This will create a localhost of a bingo card, that is set to to the Week of the year. This means that everyone that acesses this bingo card in the same week, will get the same. If you have a bingo, please take a screenshot and add it to your pull request review.