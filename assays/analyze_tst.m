function result = analyze_tst(pairs, outputDir, options)
%ANALYZE_TST Tail suspension test weighted immobility analysis.

result = analyze_weighted_immobility(pairs, outputDir, options, 'TST', 'Summary_TST.xlsx');
end
