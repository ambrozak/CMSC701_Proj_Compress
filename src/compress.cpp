#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <cstring>
#include <algorithm>

using namespace std;

// Read the first DNA sequence from a FASTA file
string loadFirstSequence(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    string line, sequence;
    bool inSequence = false;

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (inSequence) break; // Stop after first sequence
            inSequence = true;
            continue;
        }
        if (inSequence) sequence += line;
    }

    file.close();
    return sequence;
}

int main(int argc, char* argv[]) {

    if (argc != 4) {
        printf("Use: compress in_ref.fa bits_per_char out_compressed.fa\n");
        return 1;
    }

    string inFile  = argv[1];
    int bitsPerChar = stoi(argv[2]);
    string outFile = argv[3];

    // --- Load ---
    string sequence = loadFirstSequence(inFile);
    if (sequence.empty()) {
        cerr << "Error: No sequence found in " << inFile << endl;
        return 1;
    }
    cout << "Loaded sequence of length " << sequence.size() << endl;
    int starting_size = sequence.size();

    // Uppercase the sequence
    transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);

    // --- Compress ---

    // Track which characters are currently in use.
    // DNA starts with A, T, G, C (and possibly N, etc.)
    // We'll track used chars and assign new ones from a pool.
    // Max characters allowed = 2^bitsPerChar
    int maxChars = 1 << bitsPerChar;

    // Build initial set of used characters
    map<char, bool> usedChars;
    for (char c : sequence) usedChars[c] = true;

    // Rules: maps new_char -> pair of chars it replaced
    vector<tuple<char, char, char>> rules; // (new_char, left, right)

    // Pool of candidate characters to assign as new symbols.
    // Use printable ASCII that aren't already in use.
    // We'll pick from a set that's unlikely to be in DNA sequences.
    auto getNextFreeChar = [&]() -> char {
        // Prefer characters in ranges that are clearly non-DNA
        string candidates = "!#$%&'()*+,-./:;<=>?@[]^_`{|}~"
                            "0123456789"
                            "EFHIJKLMNOPQRSUVWXYZ"
                            "abcdefghijklmnopqrstuvwxyz";
        for (char c : candidates) {
            if (!usedChars[c]) return c;
        }
        cerr << "Error: Ran out of printable ASCII characters for new symbols." << endl;
        exit(1);
    };

    // RePair loop
    while (true) {
        // Count all adjacent pairs
        map<pair<char,char>, int> pairCount;
        for (size_t i = 0; i + 1 < sequence.size(); i++) {
            pairCount[{sequence[i], sequence[i+1]}]++;
            if(sequence[i] == sequence[i+1] && sequence[i+1] == sequence[i+2]) {i++;}
        }

        if (pairCount.empty()) break;

        // Find most frequent pair
        pair<char,char> bestPair;
        int bestCount = 0;
        for (auto& [p, cnt] : pairCount) {
            if (cnt > bestCount) {
                bestCount = cnt;
                bestPair = p;
            }
        }

        // Only compress if the pair appears more than once (otherwise no gain)
        if (bestCount < 2) {
            cout << "Most frequent pair appears only once — stopping early." << endl;
            break;
        }

        // Check if we have room for a new character
        if ((int)usedChars.size() >= maxChars) {
            cout << "Reached max character limit (" << maxChars
                 << ") for " << bitsPerChar << " bits — stopping." << endl;
            break;
        }

        // Assign a new character for this pair
        char newChar = getNextFreeChar();
        usedChars[newChar] = true;
        rules.emplace_back(newChar, bestPair.first, bestPair.second);

        cout << "Rule " << rules.size() << ": '"  << newChar << "' -> '"
             << bestPair.first << "' + '" << bestPair.second
             << "'  (count=" << bestCount << ", symbols used=" << usedChars.size() << "/" << maxChars << ")" << endl;

        // Replace all non-overlapping occurrences of bestPair in sequence (left to right)
        string newSeq;
        newSeq.reserve(sequence.size());
        size_t i = 0;
        while (i < sequence.size()) {
            if (i + 1 < sequence.size() &&
                sequence[i] == bestPair.first &&
                sequence[i+1] == bestPair.second) {
                newSeq += newChar;
                i += 2;
            } else {
                newSeq += sequence[i];
                i++;
            }
        }
        sequence = move(newSeq);

    }

    cout << "Starting sequence length: " << starting_size << endl;
    cout << "Compressed sequence length: " << sequence.size() << endl;
    int ending_size = sequence.size();
    cout  << "Starting sequence length (bits): " << starting_size * 2 << endl;
    cout << "Ending sequence length (bits): " << ending_size * bitsPerChar << endl;
    cout << "Total rules: " << rules.size() << endl;

    // --- Write out ---
    ofstream out(outFile);
    if (!out.is_open()) {
        cerr << "Error: Cannot open output file " << outFile << endl;
        return 1;
    }

    // Write rules section
    out << ">RULES count=" << rules.size() << " bits_per_char=" << bitsPerChar << "\n";
    for (auto& [newChar, left, right] : rules) {
        out << newChar << " " << left << " " << right << "\n";
    }

    // Write compressed sequence
    out << ">COMPRESSED\n";
    // Write in lines of 80 characters (FASTA convention)
    const int LINE_LEN = 80;
    for (size_t i = 0; i < sequence.size(); i += LINE_LEN) {
        out << sequence.substr(i, LINE_LEN) << "\n";
    }

    out.close();
    cout << "Written to " << outFile << endl;

    return 0;
}