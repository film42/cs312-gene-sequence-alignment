using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab {
  class PairWiseAlign {
    
    /// <summary>
    /// Align only 5000 characters in each sequence.
    /// </summary>
    private int MaxCharactersToAlign = 5000; 

    private List<int> resultSet1;
    private List<int> resultSet2;

    private Dictionary<int, Dictionary<int, int>> cache = new Dictionary<int, Dictionary<int, int>>();

    private int indel = 5;
    private int sub = 1;
    private int match = -3;

    /// <summary>
    /// this is the function you implement.
    /// </summary>
    /// <param name="sequenceA">the first sequence</param>
    /// <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
    /// <param name="resultTableSoFar">the table of alignment results that has been generated so far using pair-wise alignment</param>
    /// <param name="rowInTable">this particular alignment problem will occupy a cell in this row the result table.</param>
    /// <param name="columnInTable">this particular alignment will occupy a cell in this column of the result table.</param>
    /// <returns>the alignment score for sequenceA and sequenceB.  The calling function places the result in entry rowInTable,columnInTable
    /// of the ResultTable</returns>
    public int Align(GeneSequence sequenceA, GeneSequence sequenceB, ResultTable resultTableSoFar, int rowInTable, int columnInTable) {
      // Check for equal alignment
      if (rowInTable == columnInTable) return 0;
      
      // Check Cache
      if (InCache(rowInTable, columnInTable)) return cache[rowInTable][columnInTable];

      // Linear space requirement
      resultSet1 = new List<int>(); // Prev Row
      resultSet2 = new List<int>(); // Current Row

      // New Sequences
      string seqA = "0" + sequenceA.Sequence;
      string seqB = "0" + sequenceB.Sequence;

      // Clean up the length to cap at 5000
      int seqALength = seqA.Length;
      int seqBLength = seqB.Length;
      if (seqA.Length > MaxCharactersToAlign + 1) seqALength = MaxCharactersToAlign + 1;
      if (seqB.Length > MaxCharactersToAlign + 1) seqBLength = MaxCharactersToAlign + 1;

      // Core iteration
      for (int i = 0; i < seqALength; i++) {
        for (int j = 0; j < seqBLength; j++) {
          int cost = 0;
          if (i == 0 && j == 0) cost = 0;
          else if (i == 0 && j > 0) cost = resultSet2[j - 1] + indel;
          else if (i > 0 && j == 0) cost = resultSet1[j] + indel;
          else if ((i > 0 && j > 0) && (seqA[i] == seqB[j])) {
            // match or indel
            int top = resultSet1[j] + indel;
            int left = resultSet2[j - 1] + indel;
            int diag = resultSet1[j - 1] + match;
            cost = GetSmallest(top, left, diag);
          } else if ((i > 0 && j > 0) && (seqA[i] != seqB[j])) {
            // sub or indel
            int top = resultSet1[j] + indel;
            int left = resultSet2[j - 1] + indel;
            int diag = resultSet1[j - 1] + sub;
            cost = GetSmallest(top, left, diag);
          }
          resultSet2.Add(cost);
        }

        resultSet1 = resultSet2;
        resultSet2 = new List<int>();
      }

      int result = resultSet1[resultSet1.Count - 1];

      AddToCache(rowInTable, columnInTable, result);
      AddToCache(columnInTable, rowInTable, result);

      return result;     
    }

    private void AddToCache(int row, int col, int cost) {
      // If the Y key doesn't exist, make a new Dictionary and settle it
      if (!cache.ContainsKey(row))
        cache.Add(row, new Dictionary<int, int>());

      // Otherwise just settle it
      cache[row][col] = cost;
    }

    private bool InCache(int row, int col) {
      // If the key exists for both, it's found
      if (cache.ContainsKey(row))
        if (cache[row].ContainsKey(col))
          return true;

      // Otherwise, not found
      return false;
    }

    private int GetSmallest(int a, int b, int c) {
      int smallest = a;
      if (b < smallest) smallest = b;
      if (c < smallest) smallest = c;
      return smallest;
    }
  }
}
