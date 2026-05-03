#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
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

    if (argc != 3) {
        printf("Use: compress in_ref.fa out_compressed.txt\n");
        return 1;
    }

    string inFile  = argv[1];
    string outFile = argv[2];

    // --- Load ---
    string rawSequence = loadFirstSequence(inFile);
    if (rawSequence.empty()) {
        cerr << "Error: No sequence found in " << inFile << endl;
        return 1;
    }

    // Uppercase
    transform(rawSequence.begin(), rawSequence.end(), rawSequence.begin(), ::toupper);

    cout << "Loaded sequence of length " << rawSequence.size() << endl;
    size_t startingSize = rawSequence.size();

    // --- Map base characters to integers ---
    map<char, int> charToSymbol;
    map<int, char> symbolToChar; // for the alphabet section
    int nextSymbol = 0;

    vector<int> sequence;
    sequence.reserve(rawSequence.size());

    for (char c : rawSequence) {
        if (charToSymbol.find(c) == charToSymbol.end()) {
            charToSymbol[c] = nextSymbol;
            symbolToChar[nextSymbol] = c;
            nextSymbol++;
        }
        sequence.push_back(charToSymbol[c]);
    }

    cout << "Base alphabet size: " << nextSymbol << " symbols" << endl;

    // Rules: (newSymbol, left, right)
    vector<tuple<int, int, int>> rules;

    // --- RePair loop ---
    while (true) {
        // Count all adjacent pairs
        map<pair<int,int>, int> pairCount;
        for (size_t i = 0; i + 1 < sequence.size(); i++) {
            pairCount[{sequence[i], sequence[i+1]}]++;
        }

        if (pairCount.empty()) break;

        // Find most frequent pair
        pair<int,int> bestPair;
        int bestCount = 0;
        for (auto& [p, cnt] : pairCount) {
            if (cnt > bestCount) {
                bestCount = cnt;
                bestPair = p;
            }
        }

        // Stop if no pair appears more than once (no gain possible)
        if (bestCount < 2) {
            cout << "Most frequent pair appears only once — stopping." << endl;
            break;
        }

        // Assign a new integer symbol
        int newSymbol = nextSymbol++;
        rules.emplace_back(newSymbol, bestPair.first, bestPair.second);

        int totalSymbols = nextSymbol;
        int curBits = totalSymbols <= 1 ? 1 : (int)ceil(log2((double)totalSymbols));

        cout << "Rule " << rules.size() << ": " << newSymbol
             << " -> " << bestPair.first << " + " << bestPair.second
             << "  (count=" << bestCount
             << ", symbols=" << totalSymbols
             << ", bits/sym=" << curBits
             << ", seq_bits=" << (long long)sequence.size() * curBits
             << ", seq_chars=" << (long long)sequence.size()
             << ")" << endl;

        // Replace all non-overlapping occurrences left to right
        vector<int> newSeq;
        newSeq.reserve(sequence.size());
        size_t i = 0;
        while (i < sequence.size()) {
            if (i + 1 < sequence.size() &&
                sequence[i]   == bestPair.first &&
                sequence[i+1] == bestPair.second) {
                newSeq.push_back(newSymbol);
                i += 2;
            } else {
                newSeq.push_back(sequence[i]);
                i++;
            }
        }
        sequence = move(newSeq);
    }

    // --- Summary ---
    int totalSymbols = nextSymbol;
    int bitsPerSym = totalSymbols <= 1 ? 1 : (int)ceil(log2((double)totalSymbols));

    cout << "\n--- Compression Summary ---" << endl;
    cout << "Starting length (chars):      " << startingSize << endl;
    cout << "Starting size (bits, 2bpc):   " << startingSize * 2 << endl;
    cout << "Compressed length (symbols):  " << sequence.size() << endl;
    cout << "Total symbols in alphabet:    " << totalSymbols << endl;
    cout << "Bits per symbol:              " << bitsPerSym << endl;
    cout << "Compressed size (bits):       " << (long long)sequence.size() * bitsPerSym << endl;
    cout << "Total rules:                  " << rules.size() << endl;

    // --- Write output ---
    ofstream out(outFile);
    if (!out.is_open()) {
        cerr << "Error: Cannot open output file " << outFile << endl;
        return 1;
    }

    // Alphabet section: maps integer ID -> original character
    out << ">ALPHABET count=" << charToSymbol.size() << "\n";
    for (auto& [id, c] : symbolToChar) {
        out << id << " " << c << "\n";
    }

    // Rules section: newSymbol left right, all as integers
    out << ">RULES count=" << rules.size() << "\n";
    for (auto& [newSym, left, right] : rules) {
        out << newSym << " " << left << " " << right << "\n";
    }

    // Compressed sequence: space-separated integers, 20 per line for readability
    out << ">COMPRESSED length=" << sequence.size()
        << " total_symbols=" << totalSymbols
        << " bits_per_symbol=" << bitsPerSym << "\n";
    const int PER_LINE = 20;
    for (size_t i = 0; i < sequence.size(); i++) {
        out << sequence[i];
        if ((i + 1) % PER_LINE == 0 || i + 1 == sequence.size()) {
            out << "\n";
        } else {
            out << " ";
        }
    }

    out.close();
    cout << "Written to " << outFile << endl;

    return 0;
}