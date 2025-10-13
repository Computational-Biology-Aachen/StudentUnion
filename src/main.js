import './style.css';

const bingoPhrases = [
    "Code runs without errors on first try", "Excellent README.md with run instructions",
    "Dependencies are listed in a requirements file", "No hard-coded file paths",
    "Descriptive variable names (e.g., 'gene_counts' not 'd')", "Helpful comments explaining the 'why'",
    "Function has a clear docstring", "Plot has a title and labeled axes",
    "Code is organized into functions", "Avoided a slow 'for' loop on a DataFrame",
    "A 'magic number' is defined as a constant", "Learned a new library or function from this code",
    "Unused code/imports are removed", "Consistent code formatting (e.g., PEP 8)",
    "Random seed is set for reproducibility", "Output files have informative names",
    "The logic is easy to follow", "Error handling is present (e.g., try-except)",
    "A function is less than 30 lines long", "Used a built-in function instead of re-implementing it",
    "Code is reasonably efficient", "The file structure is well-organized",
    "Suggested a better way to visualize data", "Found and fixed a minor bug",
    "Documentation includes example usage", "No large data files committed to git",
    "Variable names are consistent", "The main script is easy to identify and run"
];

const boardBody = document.getElementById('bingo-body');
const resetButton = document.getElementById('reset-button');
const weekDisplayText = document.getElementById('week-display-text');

const GRID_SIZE = 5;
const MIDDLE_INDEX = Math.floor(GRID_SIZE * GRID_SIZE / 2);

// --- SEED GENERATION BASED ON WEEK ---

/**
 * Gets the ISO 8601 week number and year for a date.
 * @param {Date} d The date.
 * @returns {[number, number]} An array containing [year, weekNumber].
 */
function getWeekNumber(d) {
    d = new Date(Date.UTC(d.getFullYear(), d.getMonth(), d.getDate()));
    d.setUTCDate(d.getUTCDate() + 4 - (d.getUTCDay() || 7));
    const yearStart = new Date(Date.UTC(d.getUTCFullYear(), 0, 1));
    const weekNo = Math.ceil((((d - yearStart) / 86400000) + 1) / 7);
    return [d.getUTCFullYear(), weekNo];
}

/**
 * A simple seeded pseudorandom number generator (PRNG).
 * @param {number} seed The seed.
 * @returns {function(): number} A function that returns a random number between 0 and 1.
 */
function mulberry32(seed) {
    return function() {
      let t = seed += 0x6D2B79F5;
      t = Math.imul(t ^ t >>> 15, t | 1);
      t ^= t + Math.imul(t ^ t >>> 7, t | 61);
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    }
}

// Generate the seed for the current week
const today = new Date();
const [year, weekNumber] = getWeekNumber(today);
const weekSeed = parseInt(`${year}${weekNumber}`); // e.g., 202542
const seededRandom = mulberry32(weekSeed); // Our custom random function

// --- BOARD GENERATION ---

/**
 * Shuffles an array using a provided random function.
 * @param {Array} array The array to shuffle.
 * @param {function(): number} randomFunc The PRNG function.
 * @returns {Array} The shuffled array.
 */
function shuffle(array, randomFunc) {
    const newArray = [...array];
    for (let i = newArray.length - 1; i > 0; i--) {
        const j = Math.floor(randomFunc() * (i + 1));
        [newArray[i], newArray[j]] = [newArray[j], newArray[i]];
    }
    return newArray;
}

function generateBoard() {
    // Shuffle the master list of phrases using our seeded random function
    const shuffledPhrases = shuffle([...bingoPhrases], seededRandom);
    const selectedPhrases = shuffledPhrases.slice(0, GRID_SIZE * GRID_SIZE);
    
    // Set the free space
    selectedPhrases[MIDDLE_INDEX] = "FREE SPACE (Code Review)";

    // Render the board HTML
    boardBody.innerHTML = '';
    selectedPhrases.forEach((phrase, index) => {
        if (index % GRID_SIZE === 0) {
            boardBody.insertAdjacentHTML('beforeend', '<tr></tr>');
        }
        const row = boardBody.lastChild;
        row.insertAdjacentHTML('beforeend', `<td id="cell-${index}">${phrase}</td>`);
    });

    // Auto-mark the free space
    document.getElementById(`cell-${MIDDLE_INDEX}`).classList.add('marked');

    loadState();
}

// --- STATE MANAGEMENT ---

const STORAGE_KEY_STATE = `bingoState-${weekSeed}`; // State is unique to the week's card

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

// --- EVENT LISTENERS ---

boardBody.addEventListener('click', (event) => {
    if (event.target.tagName === 'TD') {
        event.target.classList.toggle('marked');
        saveState();
    }
});

resetButton.addEventListener('click', () => {
    localStorage.removeItem(STORAGE_KEY_STATE); // Clear saved marks
    // Visually un-mark all squares and then re-render the initial state
    const cells = document.querySelectorAll('#bingo-body td');
    cells.forEach(cell => cell.classList.remove('marked'));
    document.getElementById(`cell-${MIDDLE_INDEX}`).classList.add('marked');
    saveState();
});


// --- INITIAL LOAD ---

// Update the display text
weekDisplayText.textContent = `Weekly Card (Week ${weekNumber}, ${year})`;
// Generate the board
generateBoard();