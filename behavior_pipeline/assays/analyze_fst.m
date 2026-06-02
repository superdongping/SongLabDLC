function result = analyze_fst(pairs, outputDir, options)
%ANALYZE_FST Forced swim test weighted immobility analysis.

result = analyze_weighted_immobility(pairs, outputDir, options, 'FST', 'Summary_FST.xlsx');
end
