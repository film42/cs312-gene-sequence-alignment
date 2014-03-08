using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab {

  class GeneNode {
    public GeneNode Prev;
    public string X;
    public string Y;
    public int Cost;

    public GeneNode(GeneNode prev, char x, char y, int cost) {
      this.Prev = prev;
      this.X = x.ToString();
      this.Y = y.ToString();
      this.Cost = cost;
    }
  }


  class GeneSequenceAlign {

    private int MaxCharactersToAlign = 100;

    private int indel = 5;
    private int sub = 1;
    private int match = -3;

    public GeneNode Align(GeneSequence sequenceA, GeneSequence sequenceB, ResultTable resultTableSoFar, int rowInTable, int columnInTable) {
      // Check for equal alignment
      if (rowInTable == columnInTable) return new GeneNode(null, '0', '0', 0);

      // Linear space requirement
      var resultSet1 = new List<GeneNode>(); // Prev Row
      var resultSet2 = new List<GeneNode>(); // Current Row

      // New Sequences
      string seqA = "0" + sequenceA.Sequence;
      string seqB = "0" + sequenceB.Sequence;

      // Clean up the length to cap at 5000
      int seqALength = seqA.Length;
      int seqBLength = seqB.Length;
      if (seqA.Length > MaxCharactersToAlign + 1) seqALength = MaxCharactersToAlign + 1;
      if (seqB.Length > MaxCharactersToAlign + 1) seqBLength = MaxCharactersToAlign + 1;

      // Core alignment algorithm
      for (int i = 0; i < seqALength; i++) {
        for (int j = 0; j < seqBLength; j++) {
          GeneNode acc = null;
          // Starting position
          if (i == 0 && j == 0) acc = new GeneNode(null, '0', '0', 0);
          // Edge case
          else if (i == 0 && j > 0)
            // Get the node from the left cell, and add as indel
            acc = new GeneNode(resultSet2[j - 1], '-', seqB[j], (resultSet2[j - 1].Cost + indel));
          // Edge case
          else if (i > 0 && j == 0)
            // Get the node from the top cell, and add as indel
            acc = new GeneNode(resultSet1[j], seqA[i], '-', (resultSet1[j].Cost + indel));

          // Match Case
          else if ((i > 0 && j > 0) && (seqA[i] == seqB[j])) {
            // Match or indel
            var top = new GeneNode(resultSet1[j], seqA[i], '-', (resultSet1[j].Cost + indel));
            var left = new GeneNode(resultSet2[j - 1], '-', seqB[j], (resultSet2[j - 1].Cost + indel));
            var diag = new GeneNode(resultSet1[j - 1], seqA[i], seqB[j], (resultSet1[j - 1].Cost + match));
            // Get node neighbor with smallest cost in order of: left, top, diag
            acc = GetSmallest(left, top, diag);
          } 

          // Subsitution Case
          else if ((i > 0 && j > 0) && (seqA[i] != seqB[j])) {
            // Sub or indel
            var top = new GeneNode(resultSet1[j], seqA[i], '-', (resultSet1[j].Cost + indel));
            var left = new GeneNode(resultSet2[j - 1], '-', seqB[j], (resultSet2[j - 1].Cost + indel));
            var diag = new GeneNode(resultSet1[j - 1], seqA[i], seqB[j], (resultSet1[j - 1].Cost + sub));
            // Get node neighbor with smallest cost in order of: left, top, diag
            acc = GetSmallest(left, top, diag);
          }
          // Add to lower row (current row) results
          resultSet2.Add(acc);
        }

        // Make lower row the new upper row and clear the old lower row
        resultSet1 = resultSet2;
        resultSet2 = new List<GeneNode>();
      }
      // return the node at furthest (col, row) from origin
      return  resultSet1[resultSet1.Count - 1];
    }

    private GeneNode GetSmallest(GeneNode a, GeneNode b, GeneNode c) {
      var smallest = a;
      if (b.Cost < smallest.Cost) smallest = b;
      if (c.Cost < smallest.Cost) smallest = c;
      return smallest;
    }

  }


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
          // First cell, no cost value
          if (i == 0 && j == 0) cost = 0;
          // Get the cost from the left cell, and add as indel
          else if (i == 0 && j > 0) cost = resultSet2[j - 1] + indel;
          // Get the cost from the top cell, and add as indel
          else if (i > 0 && j == 0) cost = resultSet1[j] + indel;
          else if ((i > 0 && j > 0) && (seqA[i] == seqB[j])) {
            // Match or indel
            int top = resultSet1[j] + indel;
            int left = resultSet2[j - 1] + indel;
            int diag = resultSet1[j - 1] + match;
            // Get smallest costs in order of: top, left, diag
            cost = GetSmallest(top, left, diag);
          } else if ((i > 0 && j > 0) && (seqA[i] != seqB[j])) {
            // Sub or indel
            int top = resultSet1[j] + indel;
            int left = resultSet2[j - 1] + indel;
            int diag = resultSet1[j - 1] + sub;
            // Get smallest costs in order of: top, left, diag
            cost = GetSmallest(top, left, diag);
          }
          // Add to the lower row.. the result set.
          resultSet2.Add(cost);
        }

        // Make the lower row the new upper row
        resultSet1 = resultSet2;
        // Clear the old lower row
        resultSet2 = new List<int>();
      }

      // Get result
      int result = resultSet1[resultSet1.Count - 1];

      // Add to the cache
      AddToCache(rowInTable, columnInTable, result);
      AddToCache(columnInTable, rowInTable, result);

      // Return the new calculated cost
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
