document.addEventListener('DOMContentLoaded', () => {

    const bingoPhrases = [
        "Code runs without errors on first try",
        "Excellent README.md with run instructions",
        "Dependencies are listed in a requirements file",
        "No hard-coded file paths",
        "Descriptive variable names (e.g., 'gene_counts' not 'd')",
        "Helpful comments explaining the 'why'",
        "Function has a clear docstring",
        "Plot has a title and labeled axes",
        "Code is organized into functions",
        "Avoided a slow 'for' loop on a DataFrame",
        "A 'magic number' is defined as a constant",
        "Learned a new library or function from this code",
        "Unused code/imports are removed",
        "Consistent code formatting (e.g., PEP 8)",
        "Random seed is set for reproducibility",
        "Output files have informative names",
        "The logic is easy to follow",
        "Error handling is present (e.g., try-except)",
        "A function is less than 30 lines long",
        "Used a built-in function instead of re-implementing it",
        "Code is reasonably efficient",
        "The file structure is well-organized",
        "Suggested a better way to visualize data",
        "Found and fixed a minor bug",
        "Documentation includes example usage",
        "No large data files committed to git",
        "Variable names are consistent",
        "The main script is easy to identify and run"
    ];

    const boardBody = document.getElementById('bingo-body');
    const resetButton = document.getElementById('reset-button');
    const GRID_SIZE = 5;
    const STORAGE_KEY_SQUARES = 'bingoSquares';
    const STORAGE_KEY_STATE = 'bingoState';

    // Fisher-Yates shuffle algorithm to randomize an array
    function shuffle(array) {
        for (let i = array.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [array[i], array[j]] = [array[j], array[i]];
        }
        return array;
    }

    function generateBoard() {
        // Check if a board is already stored for this session
        let squares = JSON.parse(localStorage.getItem(STORAGE_KEY_SQUARES));

        if (!squares) {
            squares = shuffle([...bingoPhrases]).slice(0, GRID_SIZE * GRID_SIZE);
            // Add Free Space in the middle
            const middleIndex = Math.floor(squares.length / 2);
            squares[middleIndex] = "FREE SPACE (Code Review)";
            localStorage.setItem(STORAGE_KEY_SQUARES, JSON.stringify(squares));
        }

        boardBody.innerHTML = '';
        let squareIndex = 0;
        for (let i = 0; i < GRID_SIZE; i++) {
            const row = document.createElement('tr');
            for (let j = 0; j < GRID_SIZE; j++) {
                const cell = document.createElement('td');
                cell.id = `cell-${squareIndex}`;
                cell.textContent = squares[squareIndex];
                if (squareIndex === Math.floor(GRID_SIZE * GRID_SIZE / 2)) {
                    cell.classList.add('marked'); // Auto-mark the free space
                }
                row.appendChild(cell);
                squareIndex++;
            }
            boardBody.appendChild(row);
        }
        loadState();
    }

    function saveState() {
        const state = {};
        const cells = document.querySelectorAll('#bingo-body td');
        cells.forEach(cell => {
            state[cell.id] = cell.classList.contains('marked');
        });
        localStorage.setItem(STORAGE_KEY_STATE, JSON.stringify(state));
    }

    function loadState() {
        const state = JSON.parse(localStorage.getItem(STORAGE_KEY_STATE));
        if (state) {
            Object.keys(state).forEach(cellId => {
                const cell = document.getElementById(cellId);
                if (cell && state[cellId]) {
                    cell.classList.add('marked');
                }
            });
        }
    }

    // Event listener for clicking on squares (using event delegation)
    boardBody.addEventListener('click', (event) => {
        if (event.target.tagName === 'TD') {
            event.target.classList.toggle('marked');
            saveState();
        }
    });

    // Event listener for the reset button
    resetButton.addEventListener('click', () => {
        localStorage.removeItem(STORAGE_KEY_STATE);
        localStorage.removeItem(STORAGE_KEY_SQUARES);
        generateBoard(); // This will create a brand new, random board
    });

    // Initial setup
    generateBoard();
});