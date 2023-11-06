package org.monarchinitiative.exomiser.core.genome;



import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.antlr.v4.runtime.tree.Tree;
import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.jannovar.JannovarDataSourceLoader;
import org.monarchinitiative.exomiser.core.model.ChromosomalRegionIndex;
import org.monarchinitiative.exomiser.core.model.VariantAnnotation;
import org.monarchinitiative.svart.CoordinateSystem;
import org.monarchinitiative.svart.GenomicVariant;
import org.monarchinitiative.svart.Strand;


import java.io.*;
import java.lang.reflect.Array;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;


class TestJannovarAnnotationOnVCF {

    private Path transcriptFilePath = Path.of("/Users/carlo/Documents/exomiser-data/2302_hg38_transcripts_ensembl.ser");

    private final JannovarVariantAnnotator jannovarVariantAnnotator = new JannovarVariantAnnotator(GenomeAssembly.HG38, JannovarDataSourceLoader.loadJannovarData(transcriptFilePath), ChromosomalRegionIndex.empty());

    private List<VariantAnnotation> annotate(VariantAnnotator instance, String contig, int start, String ref, String alt) {
        GenomicVariant variant = variant(contig, start, ref, alt);
        return instance.annotate(variant);
    }

    private GenomicVariant variant(String contig, int start, String ref, String alt) {
//        VariantTrimmer.VariantPosition variantPosition = VariantTrimmer.leftShiftingTrimmer(VariantTrimmer.retainingCommonBase()).trim(Strand.POSITIVE, start, ref, alt);
        return GenomicVariant.of(GenomeAssembly.HG38.getContigByName(contig), Strand.POSITIVE, CoordinateSystem.ONE_BASED, start, ref, alt);
    }

    @Test
    void getFailure() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch37.vcf.gz");
        Path fileAPath = Path.of("/Users/carlo/Desktop/jannovar/jannovar-cli-0.38/jannovarFromCsvMadeByClinVarVcfAnnotated_grch37.txt");

        Map<String, List<String>> fileAMap = new HashMap<>();

        try (BufferedReader bufferedReader = new BufferedReader(new FileReader(fileAPath.toString()))) {
            String lineA;
            while ((lineA = bufferedReader.readLine()) != null) {
                String[] tokensA = lineA.split(",");

                String variantA = String.format("%s-%s-%s-%s", tokensA[0].split("chr")[1], tokensA[1], tokensA[2], tokensA[3]);
                String geneSymbolA = tokensA[4].split(":")[0];
                String variantEffectA = tokensA[5];
                List<String> arrayList = List.of(geneSymbolA, variantEffectA);
                fileAMap.put(variantA, arrayList);

            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        Set<String> uniqueClinVarSet = new TreeSet<>();
        Set<String> allClinVarSet = new TreeSet<>();
        Set<String> partialJannoSet = new TreeSet<>();
        Set<String> allJannoSet = new TreeSet<>();
        TreeMap<String, Integer> uniqueEffectCounts = new TreeMap<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/delete--testRealJannovar.txt"));
        BufferedWriter effectWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/delete--uniqueEffectsJanno.txt"));
        BufferedWriter uniqueClinVarWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/delete--partialSetClinVarMc.txt"));
        BufferedWriter allClinVarMcWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/delete--allSetClinVarMc.txt")); // whole file
        BufferedWriter partJannoMcWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/delete--partialSetJannoMc.txt")); // where match = 0 (but all, no mc.contains)
        BufferedWriter allJannoMcWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/delete--allSetJannoMc_2.txt")); // whole file
        BufferedWriter catchMissenseSynonymous = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/delete__missenseSynonymous.txt"));
        BufferedWriter lmf1writer = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/__LMF1.txt")))
        {
            int totalVariants = 0;
            int totalVariantsWithMultipleConsequences = 0;
            int totalMatchesFromVariantsWithMultipleConsequences = 0;
            String line;
            Map<String, Map<String, Integer>> lmf1Stats = new HashMap<>();

            int isMatch = 0;
            Map<String, String> infoMap = new HashMap<>();
            bw.write("Variant\tclinvarGeneSymbol\tclinvarFirstEffect\tjannovarGeneSymbol\tjannovarFirstEffect\tisMatch\n");
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                totalVariants++;
                String[] tokens = line.split("\t");
                String[] infoTokens = tokens[7].split(";");
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }
                String chrom = tokens[0];
                int pos = Integer.parseInt(tokens[1]);
                String ref = tokens[3];
                String alt = tokens[4];
                String mc = infoMap.get("MC");
                String geneInfo = infoMap.get("GENEINFO");


                if (mc != null && geneInfo != null) {
                    List<VariantAnnotation> variantAnnotations = annotate(jannovarVariantAnnotator, chrom, pos, ref, alt);
                    totalVariantsWithMultipleConsequences++;

                    String clinvarGeneSymbol = geneInfo.split(":")[0];
                    String clinvarFirstEffect = mc.split(",")[0].split("\\|")[1].toLowerCase();
//                    String [] clinVarEffects = mc.split(",");
                    String[] clinVarEffects = Arrays.stream(mc.split(","))
                            .map(s -> s.replaceAll(".*\\|", "").trim())
                            .toArray(String[]::new);
                    int numberClinVarEffects = clinVarEffects.length;
                    allClinVarSet.add(clinvarFirstEffect);
                    String variant = String.format("%s-%d-%s-%s", chrom, pos, ref, alt);
                    if ( !variantAnnotations.isEmpty()) {
                        List<VariantEffect> jannVariantEffects = variantAnnotations.stream()
                                .map(va -> va.getVariantEffect())
                                .collect(Collectors.toList());
                        VariantAnnotation variantAnnotation = variantAnnotations.get(0);
                        String jannVariantEffect = variantAnnotation.getVariantEffect().getSequenceOntologyTerm();
                        allJannoSet.add(jannVariantEffect);
                        String jannGeneSymbol = variantAnnotation.getGeneSymbol();

                        if (clinvarGeneSymbol.equals("LMF1")){


                            String clinicalSignificance = infoMap.get("CLNSIG");
                            System.out.println(clinicalSignificance);

                            Map<String, Integer> clinSigCountMap = lmf1Stats.getOrDefault(clinvarFirstEffect, new HashMap<>());
                            clinSigCountMap.put(clinicalSignificance, clinSigCountMap.getOrDefault(clinicalSignificance, 0) + 1);
                            lmf1Stats.put(clinvarFirstEffect, clinSigCountMap);



                        }

                        if (jannVariantEffect.equalsIgnoreCase(clinvarFirstEffect)){
                            isMatch = 1;
                            totalMatchesFromVariantsWithMultipleConsequences++;
                        } else {
                            isMatch = 0;
                            String transformedClinvarFirsteffectString = transformClinvarFirstEffect(clinvarFirstEffect);
                            uniqueEffectCounts.put(transformedClinvarFirsteffectString + "-" + jannVariantEffect, uniqueEffectCounts.getOrDefault(transformedClinvarFirsteffectString + "-" + jannVariantEffect, 0)+1);
                            uniqueClinVarSet.add(transformedClinvarFirsteffectString);
                            partialJannoSet.add(jannVariantEffect);
                        }
                        bw.write(String.format("%s\t%s\t%s\t%s\t%s\t%d\n", variant, clinvarGeneSymbol, clinvarFirstEffect, jannGeneSymbol, jannVariantEffect, isMatch));
                        if (numberClinVarEffects > 1 && isMatch == 0 && (clinvarFirstEffect.contains("missense") && jannVariantEffect.contains("synonymous"))) {
                            catchMissenseSynonymous.write("Variant\tclinvarGeneSymbol\tclinvarFirstEffects\tjannovarGeneSymbol\tjannovarFirstEffect\n");
                            catchMissenseSynonymous.write(String.format("%s\t%s\t%s\t%s\t%s\n", variant, clinvarGeneSymbol, Arrays.toString(clinVarEffects), jannGeneSymbol, jannVariantEffects.toString()));
                        }
                    }else
                    bw.write(String.format("%s\t%s\t%s\tNA\tNA\tNA\n", variant, clinvarGeneSymbol, clinvarFirstEffect));
                }
            }

            lmf1writer.write("Gene Table for LMF1:\n");

            writeLMF1Stats(lmf1Stats, lmf1writer);


            System.out.println("totalVariantsWithMultipleConsequences :  " + totalVariantsWithMultipleConsequences);
            System.out.println("totalMatchesFromVariantsWithMultipleConsequences : " + totalMatchesFromVariantsWithMultipleConsequences);

            double totalMatchesx100 = totalMatchesFromVariantsWithMultipleConsequences * 100;
            double percentageMatches = totalMatchesx100 / totalVariantsWithMultipleConsequences;

            System.out.println("Percentage of isMatch variants: " + percentageMatches + "%");

            System.out.println("totalVariants :  " + totalVariants);
            System.out.println("totalVariantsWithMultipleConsequences :  " + totalVariantsWithMultipleConsequences);

            double totalMatchesFromVariantsWithMultipleConsx100 = totalVariantsWithMultipleConsequences * 100;
            double multipleConsVariantsToTotalVariants = totalMatchesFromVariantsWithMultipleConsx100 / totalVariants;

            System.out.println("total Percentage of multipleConsVariantsToTotalVariants : " + multipleConsVariantsToTotalVariants + "%");

            List<Map.Entry<String, Integer>> sortedEntries = new ArrayList<>(uniqueEffectCounts.entrySet());

            sortedEntries.sort((e1, e2) -> {
                String[] parts1 = e1.getKey().split("-");
                String[] parts2 = e2.getKey().split("-");

                // Compare the prefixes
                int cmp = parts1[0].compareTo(parts2[0]);

                // If prefixes are different, sort by prefix
                if (cmp != 0) {
                    return cmp;
                }
                // If prefixes are the same, sort by count in descending order
                return e2.getValue().compareTo(e1.getValue());
            });


            for (Map.Entry<String, Integer> entry : sortedEntries) {
                effectWriter.write(entry.getKey() + " Count:" + entry.getValue());
                effectWriter.newLine();
                }

            for (String s : uniqueClinVarSet){
                uniqueClinVarWriter.write(s);
                uniqueClinVarWriter.newLine();
            }

            for (String y : allClinVarSet) {
                allClinVarMcWriter.write(y);
                allClinVarMcWriter.newLine();
            }

            for (String l : partialJannoSet) {
                partJannoMcWriter.write(l);
                partJannoMcWriter.newLine();
            }

            for (String n: allJannoSet) {
                allJannoMcWriter.write(n);
                allJannoMcWriter.newLine();
            }


            }catch (IOException e){
            e.printStackTrace();
        }
    }

    @Test
    void putIntoDatabase() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch38.vcf.gz");
        Map<String, Map<String, Map<String, Integer>>> geneStats = new TreeMap<>();
        int totalClassifiedVariants = 0;


        Map<String, String> infoMap = new HashMap<>();

        int totalVariantsForGene = 0;
        try {
            Connection conn = initializeDatabase();
            createDatabaseTable(conn);

            try {BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
                String line;

                while ((line = br.readLine()) != null){
                    if (line.startsWith("#")) continue;
                    String[] tokens = line.split("\t");
                    String[] infoTokens = tokens[7].split(";");
                    for (String t : infoTokens) {
                        String[] kv = t.split("=");
                        if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                    }
                }

            }catch (IOException e){
                e.printStackTrace();
            }
        }catch (SQLException e) {
            e.printStackTrace();
        }

    }
    private Connection initializeDatabase() throws SQLException {
        return DriverManager.getConnection("jdbc:h2:~/gene-data-exomiser-db", "carlo_sa", "exomiser");
        }


    private void createDatabaseTable(Connection conn) throws SQLException{
        Statement statement = conn.createStatement();
        statement.execute("CREATE TABLE IF NOT EXISTS gene_data_exomiser (gene VARCHAR, molecular_consequence VARCHAR, likely_benign INT, benign INT, uncertain_significance INT, likely_pathogenic INT, pathogenic INT)");

    }

    private void insertGeneData(Connection conn, String gene, String molecular_consequence, int likely_benign, int benign, int uncertain_significance, int likely_pathogenic, int pathogenic ) throws SQLException {
        PreparedStatement stmt = conn.prepareStatement("INSERT INTO gene_data_exomiser (gene, molecular_consequence, likely_benign, benign, uncertain_significance, likely_pathogenic, pathogenic) VALUES (?, ?, ?, ?, ?, ?, ?)");
        stmt.setString(1, gene);
        stmt.setString(2, molecular_consequence);
        stmt.setInt(3, likely_benign);
        stmt.setInt(4, benign);
        stmt.setInt(5, uncertain_significance);
        stmt.setInt(6, likely_pathogenic);
        stmt.setInt(7, pathogenic);
        stmt.executeUpdate();
    }
    @Test
    void handleLMF1Both() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch38.vcf.gz");
        //Map<effect, Map<clinSig, count>>
        Map<String, Map<String, Integer>> lmf1StatsClinvar = new TreeMap<>();
        //Map<GeneSymbolMap<effect, Map<clinSig, count>>>
        Map<String, Map<String, Map<String, Integer>>> geneStats = new TreeMap<>();
        Map<String, Map<String, Integer>> lmf1StatsJannovar = new TreeMap<>();
        Map<String, Integer> variantStatsClinvar = new TreeMap<>();
        Map<String, Integer> variantStatsJannovar = new TreeMap<>();
        int totalClassifiedVariants = 0;


        Map<String, String> infoMap = new HashMap<>();

        int totalVariantsForGene = 0;
        try {
            Connection conn = initializeDatabase();
            createDatabaseTable(conn);


        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter lmf1writer = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/30102023fgfr2.txt"));
             BufferedWriter tableWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/30102023fgfr2.txt"))) {

            int counter = 0;
            String line;
            lmf1writer.write("Gene Table for LMF1:\n");
            while ((line = br.readLine()) != null && counter < 100) {
                counter ++;
                if (line.startsWith("#")) continue;
                String[] tokens = line.split("\t");
                String[] infoTokens = tokens[7].split(";");
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }
                //jann
                String chrom = tokens[0];
                int pos = Integer.parseInt(tokens[1]);
                String ref = tokens[3];
                String alt = tokens[4];
                //clin
                String mc = infoMap.get("MC");
                String geneInfo = infoMap.get("GENEINFO");

                if (mc != null && geneInfo != null) {

                    // clinvar process
                    String clinvarGeneSymbol = geneInfo.split(":")[0];
                    String clinvarFirstEffect = mc.split(",")[0].split("\\|")[1].toLowerCase();

//                        if (clinvarGeneSymbol.equals("RBM8A")) {
                        totalVariantsForGene++;
                        //jannovar process
                        List<VariantAnnotation> variantAnnotations = annotate(jannovarVariantAnnotator, chrom, pos, ref, alt);

                        VariantAnnotation variantAnnotation = variantAnnotations.get(0);
                        String jannVariantEffect = variantAnnotation.getVariantEffect().getSequenceOntologyTerm();
                        String jannGeneSymbol = variantAnnotation.getGeneSymbol();
                        String clinicalSignificance = infoMap.get("CLNSIG");
                        // "clinSigCountMap1:
                        Map<String, Integer> clinSigCountMap = lmf1StatsClinvar.getOrDefault(clinvarFirstEffect, new HashMap<>());
                        System.out.println("clinSigCountMap1: " + clinSigCountMap.toString());

                        Map<String, Integer> clinSigCountMapJannovar = lmf1StatsJannovar.getOrDefault(jannVariantEffect, new HashMap<>());
                        System.out.println("clinSigCountMapJannovar1: " + clinSigCountMapJannovar.toString());

                        // "clinSigCountMap2:
                        clinSigCountMap.put(clinicalSignificance, clinSigCountMap.getOrDefault(clinicalSignificance, 0) + 1);
                        System.out.println("clinSigCountMap2: " + clinSigCountMap.toString());

                        clinSigCountMapJannovar.put(clinicalSignificance, clinSigCountMapJannovar.getOrDefault(clinicalSignificance, 0) +1);
                        System.out.println("clinSigCountMapJannovar2: " + clinSigCountMapJannovar.toString());

                        lmf1StatsClinvar.put(clinvarGeneSymbol, clinSigCountMap);
                        System.out.println("lmf1StatsClinvar: " + lmf1StatsClinvar.toString());
                        lmf1StatsJannovar.put(jannVariantEffect, clinSigCountMapJannovar);
                        System.out.println("lmf1StatsJannovar: " + lmf1StatsJannovar.toString());


                        totalClassifiedVariants++;
                        variantStatsJannovar.put(clinicalSignificance, variantStatsJannovar.getOrDefault(clinicalSignificance, 0) + 1);
                        variantStatsClinvar.put(clinicalSignificance, variantStatsClinvar.getOrDefault(clinicalSignificance, 0) + 1);
                            System.out.println("lmf1StatsClinvar7: " + variantStatsClinvar.toString());
                            System.out.println("lmf1StatsJannovar7: " + variantStatsJannovar.toString());



                            for (Map.Entry<String, Map<String, Integer>> entry : lmf1StatsClinvar.entrySet()) {
                            String gene = entry.getKey();
                            Map<String, Integer> stats = entry.getValue();
                            int likely_benign = stats.getOrDefault("Likely_benign", 0);
                            int benign = stats.getOrDefault("Benign", 0);
                            int uncertain_significance = stats.getOrDefault("Uncertain_significance", 0);
                            int likely_pathogenic = stats.getOrDefault("Likely_pathogenic", 0);
                            int pathogenic = stats.getOrDefault("Pathogenic", 0);
                            System.out.println(gene + " " + clinvarFirstEffect + " " + likely_benign + " " + benign + " " + uncertain_significance + " " + likely_pathogenic + " " + pathogenic);
//                            insertGeneData(conn, gene, clinvarFirstEffect, likely_benign, benign, uncertain_significance, likely_pathogenic, pathogenic);
                        }
//                    }

                }

            }

            //clinvarstats
            writeVariantStats(lmf1writer, tableWriter, variantStatsClinvar, lmf1StatsClinvar, totalClassifiedVariants, totalVariantsForGene);
            // jannovarstats
            writeVariantStats(lmf1writer, tableWriter, variantStatsJannovar, lmf1StatsJannovar, totalClassifiedVariants, totalVariantsForGene);

        } catch (IOException e) {
            e.printStackTrace();
        }
            conn.close();

        } catch (SQLException e) {
            e.printStackTrace();
        }
    }



    private void writeVariantStats(BufferedWriter lmf1writer, BufferedWriter tableWriter, Map<String, Integer> variantStats, Map<String, Map<String, Integer>> lmf1Stats, int totalClassifiedVariants, int totalVariantsForGene) throws IOException {
        lmf1writer.write("Total classified variants: " + totalClassifiedVariants + "\n");
        int pathogenicCount = variantStats.getOrDefault("Pathogenic", 0) + variantStats.getOrDefault("Likely_pathogenic", 0);
        int benignCount = variantStats.getOrDefault("Benign", 0) + variantStats.getOrDefault("Likely_benign", 0);
        int uncertainCount = variantStats.getOrDefault("Uncertain_significance", 0);
        lmf1writer.write("Pathogenic: " + pathogenicCount + "\n");
        lmf1writer.write("Uncertain significance: " + uncertainCount + "\n");
        lmf1writer.write("Benign: " + benignCount + "\n");
        lmf1writer.newLine();
        writeStatsPercentage(lmf1Stats, lmf1writer, totalVariantsForGene);
        lmf1writer.newLine();
        writeStats(lmf1Stats, lmf1writer);
        lmf1writer.newLine();
        writeStatsTable(lmf1Stats, tableWriter, totalClassifiedVariants, variantStats);
        tableWriter.newLine();
        lmf1writer.newLine();
    }


    private void writeStatsPercentage(Map<String, Map<String, Integer>> lmf1Stats, BufferedWriter lmf1writer, int totalVariantsForGene) throws IOException {
        for (Map.Entry<String, Map<String, Integer>> entry : lmf1Stats.entrySet()) {
            lmf1writer.write("Effect: " + entry.getKey() + "\n");
            for (Map.Entry<String, Integer> subEntry : entry.getValue().entrySet()) {
                double percentage = ((double) subEntry.getValue() / totalVariantsForGene) * 100;
                lmf1writer.write("\t" + subEntry.getKey() + ": " + String.format("%.2f", percentage) + "%\n");
            }
        }
    }

    void writeStatsTable(Map<String, Map<String, Integer>> statsMap, BufferedWriter writer, int totalClassifiedVariants, Map<String, Integer> variantStatsJannoOrCLinvar) throws IOException {
        writer.write("Total classified variants: " + totalClassifiedVariants + "\n");
        int pathogenicCount = variantStatsJannoOrCLinvar.getOrDefault("Pathogenic", 0) + variantStatsJannoOrCLinvar.getOrDefault("Likely_pathogenic", 0);
        int benignCount = variantStatsJannoOrCLinvar.getOrDefault("Benign", 0) + variantStatsJannoOrCLinvar.getOrDefault("Likely_benign", 0);
        int uncertainCount = variantStatsJannoOrCLinvar.getOrDefault("Uncertain_significance", 0);
        writer.write("Pathogenic: " + pathogenicCount + "\n");
        writer.write("Uncertain significance: " + uncertainCount + "\n");
        writer.write("Benign: " + benignCount + "\n");
        writer.newLine();
        writer.write("Effect\tBenign\tUncertain_Significance\tConflicting_Interpretations_of_Pathogenicity\tBenign/Likely_Benign\tLikely_Benign\tLikely_Pathogenic\tPathogenic\tPathogenic/Likely\tPathogenic\tTotal\n");
        for (Map.Entry<String, Map<String, Integer>> entry : statsMap.entrySet()) {
            String effect = entry.getKey();
            Map<String, Integer> stats = entry.getValue();
            int total = stats.values().stream().mapToInt(Integer::intValue).sum();

            writer.write(effect);
            writer.write("\t" + stats.getOrDefault("Benign", 0) + "(" + getPercentage(stats.getOrDefault("Benign", 0), total) + "%)");
            writer.write("\t" + stats.getOrDefault("Uncertain_significance", 0)+ "(" + getPercentage(stats.getOrDefault("Uncertain_significance", 0), total) + "%)");
            writer.write("\t" + stats.getOrDefault("Conflicting_interpretations_of_pathogenicity", 0) + "(" + getPercentage(stats.getOrDefault("Conflicting_interpretations_of_pathogenicity", 0), total) + "%)");
            writer.write("\t" + stats.getOrDefault("Benign/Likely_benign", 0) + "(" + getPercentage(stats.getOrDefault("Benign/Likely_benign", 0), total) + "%)");
            writer.write("\t" + stats.getOrDefault("Likely_benign", 0)+ "(" + getPercentage(stats.getOrDefault("Likely_benign", 0), total) + "%)");
            writer.write("\t" + stats.getOrDefault("Likely_pathogenic", 0)+ "(" + getPercentage(stats.getOrDefault("Likely_pathogenic", 0), total) + "%)");
            writer.write("\t" + stats.getOrDefault("Pathogenic", 0)+ "(" + getPercentage(stats.getOrDefault("Pathogenic", 0), total) + "%)");
            writer.write("\t" + stats.getOrDefault("Likely_pathogenic/Pathogenic", 0)+ "(" + getPercentage(stats.getOrDefault("Likely_pathogenic/Pathogenic", 0), total) + "%)");
            writer.write("\t" + total );
            writer.newLine();
        }
    }

    private int getPercentage(int value, int total){
        return (total == 0) ? 0: (value * 100) / total;
    }

    @Test
    void handleLMF1() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch37.vcf.gz");
        Map<String, Map<String, Integer>> lmf1Stats = new HashMap<>();
        Map<String, Integer> variantStats = new HashMap<>();
        int totalClassifiedVariants = 0;

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter lmf1writer = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableMid/__LMF1.txt"))) {

            String line;
            lmf1writer.write("Gene Table for LMF1:\n");
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] tokens = line.split("\t");
                Map<String, String> infoMap = new HashMap<>();
                String[] infoTokens = tokens[7].split(";");
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }

                String mc = infoMap.get("MC");
                String geneInfo = infoMap.get("GENEINFO");
                if (mc != null && geneInfo != null) {
                    String clinvarGeneSymbol = geneInfo.split(":")[0];
                    String clinvarFirstEffect = mc.split(",")[0].split("\\|")[1].toLowerCase();

                    if (clinvarGeneSymbol.equals("LMF1")) {
                        String clinicalSignificance = infoMap.get("CLNSIG");
                        Map<String, Integer> clinSigCountMap = lmf1Stats.getOrDefault(clinvarFirstEffect, new HashMap<>());
                        clinSigCountMap.put(clinicalSignificance, clinSigCountMap.getOrDefault(clinicalSignificance, 0) + 1);
                        lmf1Stats.put(clinvarFirstEffect, clinSigCountMap);

                        totalClassifiedVariants++;
                        variantStats.put(clinicalSignificance, variantStats.getOrDefault(clinicalSignificance, 0) + 1);
                    }
                }
            }

            lmf1writer.write("Total classified variants: " + totalClassifiedVariants + "\n");
            int pathogenicCount = variantStats.getOrDefault("Pathogenic", 0) + variantStats.getOrDefault("Likely_pathogenic", 0);
            int benignCount = variantStats.getOrDefault("Benign", 0) + variantStats.getOrDefault("Likely_benign", 0);
            int uncertainCount = variantStats.getOrDefault("Uncertain_significance", 0);

            lmf1writer.write("Pathogenic: " + pathogenicCount + "\n");
            lmf1writer.write("Uncertain significance: " + uncertainCount + "\n");
            lmf1writer.write("Benign: " + benignCount + "\n");

            writeStats(lmf1Stats, lmf1writer);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    void writeStats(Map<String, Map<String, Integer>> lmf1Stats, BufferedWriter lmf1writer) throws IOException {
        for (Map.Entry<String, Map<String, Integer>> entry : lmf1Stats.entrySet()) {
            lmf1writer.write("Effect: " + entry.getKey() + "\n");
            for (Map.Entry<String, Integer> subEntry : entry.getValue().entrySet()) {
                lmf1writer.write("\t" + subEntry.getKey() + ": " + subEntry.getValue() + "\n");
            }
        }
    }


    private String transformClinvarFirstEffect(String clinvarFirstEffect) {
        switch (clinvarFirstEffect) {
            case "3_prime_utr_variant":
            case "5_prime_utr_variant":
                return "Regulatory and/or UTR variants";
            case "frameshift_variant":
            case "nonsense":
            case "stop_lost":
                return "Truncating Variants";
            case "missense_variant":
                return "Missense Variants";
            case "splice_acceptor_variant":
            case "splice_donor_variant":
                return "Splice Variants";
            case "non-coding_transcript_variant":
            case "intron_variant":
                return "Non-coding transcript variant";
            case "initiator_codon_variant":
            case "inframe_deletion":
            case "inframe_indel":
            case "inframe_insertion":
            case "genic_downstream_transcript_variant":
                return "Non-Truncating Coding Sequence Variants";
            case "synonymous_variant":
                return "Synonymous Variants";
            default:
                return clinvarFirstEffect;
        }
    }

    private void writeLMF1Stats(Map<String, Map<String, Integer>> lmf1Stats, BufferedWriter lmf1writer) {
        try {
            lmf1writer.write("Gene Table for LMF1:\n");
            for (Map.Entry<String, Map<String, Integer>> entry : lmf1Stats.entrySet()) {
                String variantType = entry.getKey();
                Map<String, Integer> clinicalSignificanceCounts = entry.getValue();
                for (Map.Entry<String, Integer> subEntry : clinicalSignificanceCounts.entrySet()) {
                    String clinicalSignificance = subEntry.getKey();
                    Integer count = subEntry.getValue();
                    lmf1writer.write(String.format("%s = %s, Count %d\n", variantType, clinicalSignificance, count));
                }
            }
        } catch (IOException e) {
            // Handle the exception as appropriate for your application
            e.printStackTrace();
        }
    }





    @Test
    void concatenateAndReadBothFiles() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch37.vcf.gz");
        Path fileAPath = Path.of("/Users/carlo/Desktop/jannovar/jannovar-cli-0.38/jannovarFromCsvMadeByClinVarVcfAnnotated_grch37.txt");
        Map<String, String> fileAMap = new HashMap<>();

        // transform fileA
        try (BufferedReader brA = new BufferedReader(new FileReader(fileAPath.toString()))) {
            String lineA;
            int counter = 0;
            while ((lineA = brA.readLine()) != null && counter < 1000 ) {
                counter++;
                String[] tokensA = lineA.split(",");
                String variantA = String.format("%s-%s-%s-%s", tokensA[0].split("chr")[1], tokensA[1], tokensA[2], tokensA[3]);
                System.out.println(variantA);
                String geneSymbolA = tokensA[4].split(":")[0];
                String variantEffectA = tokensA[5];
                fileAMap.put(variantA, geneSymbolA + "\t" + variantEffectA);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/testComparison.txt"));
             BufferedWriter effectWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/uniqueEffects3.txt"))) {

            bw.write("Variant\tclinvarGeneSymbol\tclinvarFirstEffect\tjannovarGeneSymbol\tjannovarFirstEffect\tisMatch\n");
            String line;
            int count = 0;
            while ((line = br.readLine()) != null && count < 5000) {
                count++;
                if (line.startsWith("#")) continue;
                String[] tokens = line.split("\t"), infoTokens = tokens[7].split(";");
                Map<String, String> infoMap = new HashMap<>();
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }
                String chrom = tokens[0];
                int pos = Integer.parseInt(tokens[1]);
                String ref = tokens[3];
                String alt = tokens[4];
                String mc = infoMap.get("MC");
                String geneInfo = infoMap.get("GENEINFO");
                if (mc != null && geneInfo != null && mc.contains(",")) {
                    String clinvarGeneSymbol = geneInfo.split(":")[0];
                    String clinvarFirstEffect = mc.split(",")[0].split("\\|")[1].toLowerCase();
                    String variant = String.format("%s-%d-%s-%s", chrom, pos, ref, alt);
                    String [] fileAInfo = fileAMap.get(variant).split("\t");

                    bw.write(String.format("%s\t%s", variant, fileAInfo.toString()));

                    if (fileAMap.containsKey(variant)) {
//                        String [] fileAInfo = fileAMap.get(variant).split("\t");
//                        System.out.println(fileAInfo);
//                        bw.write(String.format("%s\t%s\t%s\t%s\n", variant, clinvarGeneSymbol, clinvarFirstEffect, fileAInfo[0], fileAInfo[1]));
                    } else {
                        bw.write(String.format("%s\t%s\t%s\tNA\tNA\n", variant, clinvarGeneSymbol, clinvarFirstEffect));
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void annotateVariantsAndCompareNotWorkingCauseJannovarVariantAnnotatorHasNotAllTranscriptModels() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch37.vcf.gz");
        Set<String> uniqueEffects = new HashSet<>();

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/variantsComparisonGrch37.txt"));
             BufferedWriter effectWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/uniqueEffects3.txt"))) {

            bw.write("Variant\tclinvarGeneSymbol\tclinvarFirstEffect\tjannovarGeneSymbol\tjannovarFirstEffect\tisMatch\n");
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] tokens = line.split("\t"), infoTokens = tokens[7].split(";");
                Map<String, String> infoMap = new HashMap<>();
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }
                String chrom = tokens[0];
                int pos = Integer.parseInt(tokens[1]);
                String ref = tokens[3];
                String alt = tokens[4];
                String mc = infoMap.get("MC");
                String geneInfo = infoMap.get("GENEINFO");
                if (mc != null && geneInfo != null && mc.contains(",")) {
                    String clinvarGeneSymbol = geneInfo.split(":")[0];
                    String clinvarFirstEffect = mc.split(",")[0].split("\\|")[1].toLowerCase();
                    uniqueEffects.add(clinvarFirstEffect);
//                    List<VariantAnnotation> annotatedVariants = annotate(instance, chrom, pos, ref, alt);
//                    String jannovarGeneSymbol = annotatedVariants.get(0).getGeneSymbol();
//                    if (!annotatedVariants.isEmpty()) {
//                            String jannovarFirstEffect = annotatedVariants.get(0).getVariantEffect().toString().toLowerCase();
//                            int isMatch = 0;
//                            if (clinvarFirstEffect.equals(jannovarFirstEffect) || (clinvarFirstEffect.contains("intron_variant") && jannovarFirstEffect.contains("intron_variant"))) {
//                                isMatch = 1;
//                            }
                            bw.write(String.format("%s-%d-%s-%s\t%s\t%s\t%s\t%s\t%d\n", chrom, pos, ref, alt, clinvarGeneSymbol, clinvarFirstEffect)); // jannovarGeneSymbol, jannovarFirstEffect, isMatch
//                        }
                    }
                }
                for (String effect : uniqueEffects) {
                    effectWriter.write(effect);
                    effectWriter.newLine();
                }
            } catch(IOException e){
                e.printStackTrace();
            }
        }

    @Test
    void toAnnotateAcsvWithJannovarOnlyMultipleConsequenceVariants() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch37.vcf.gz");
        String outputPath = "/Users/carlo/Desktop/GeneTableStart/csvToPutIntoJannovarForAnnotationOnlyMultipleConsequences.txt";

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath))){

//            bw.write("chr,pos,ref,alt\n");
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;

                String[] tokens = line.split("\t");
                String [] infoTokens = tokens[7].split(";");
                Map<String, String> infoMap = new HashMap<>();
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }
                String token = tokens[0];
                String chr = "chr" + token;

                if (chr.contains("_")) continue;
                if (chr.contains("NT")) continue;

                int pos = Integer.parseInt(tokens[1]);
                String ref = tokens[3];
                String alt = tokens[4];
//                if (alt.contains(".")) alt = "";
                if (".".equals(alt)){
                    alt = "";
                }


                String mc = infoMap.get("MC");
                if (mc != null && mc.contains(",")) {
                    String clinvarFirstEffect = mc.split(",")[0].split("\\|")[1].toLowerCase();
                    bw.write(String.format("%s,%d,%s,%s\n", chr, pos, ref, alt));
                }
            }
            } catch(IOException e){
                e.printStackTrace();
        }
    }

    @Test
    void toAnnotateAcsvWithJannovarAllVariants() { // chr, pos, ref, alt
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvarGrch37.vcf.gz");
        String outputPath = "/Users/carlo/Desktop/GeneTableStart/csvToPutIntoJannovarForAnnotation.txt";

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter bw = new BufferedWriter(new FileWriter(outputPath))) {

            bw.write("Chr,Pos,ref,alt\n");
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;

                String[] tokens = line.split("\t");

                String chrom = tokens[0];
                int pos = Integer.parseInt(tokens[1]);
                String ref = tokens[3];
                String alt = tokens[4];
                bw.write(String.format("%s,%d,%s,%s\n", chrom, pos, ref, alt));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



    @Test
    void annotateVariantsAndCompare2() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Path outputPath = Path.of("/Users/carlo/Desktop/GeneTableStart/variantsComparison2.tsv");
        Path uniqueEffectsPath = Path.of("/Users/carlo/Desktop/GeneTableStart/uniqueEffects2.txt");
        Path annotatorOutputPath = Path.of("/Users/carlo/Desktop/GeneTableStart/annotatorOutput2.txt");

        Set<String> uniqueEffects = new HashSet<>();

        int counter = 0;

        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))));
             BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/variantsComparison2.tsv"));
             BufferedWriter effectWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/uniqueEffects2.txt"));
             BufferedWriter annotatorWriter = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/annotatorOutput2.txt"))) {

            bw.write("Variant\tclinvarGeneSymbol\tclinvarFirstEffect\tjannovarGeneSymbol\tjannovarFirstEffect\tisMatch\n");

            String line;
            while ((line = br.readLine()) != null && counter < 100) {
                if (line.startsWith("#")) continue;

                counter++;

                String[] tokens = line.split("\t"), infoTokens = tokens[7].split(";");
                Map<String, String> infoMap = new HashMap<>();
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }

                String chrom = tokens[0];
                int pos = Integer.parseInt(tokens[1]);
                String ref = tokens[3];
                String alt = tokens[4];

                String mc = infoMap.get("MC");
                String geneInfo = infoMap.get("GENEINFO");
                if (mc != null && geneInfo != null && mc.contains(",")) {
                    String clinvarGeneSymbol = geneInfo.split(":")[0];
                    String clinvarFirstEffect = mc.split(",")[0].split("\\|")[1].toLowerCase();

                    uniqueEffects.add(clinvarFirstEffect);

                    List<VariantAnnotation> annotatedVariants = annotate(jannovarVariantAnnotator, chrom, pos, ref, alt);
                    var x = annotatedVariants.get(0).getTranscriptAnnotations().toString();
                    String jannovarGeneSymbol = annotatedVariants.get(0).getGeneSymbol();
                    String jannovarFirstEffect = annotatedVariants.get(0).getVariantEffect().toString().toLowerCase();

                    // Write the annotation results to the output file
                    annotatorWriter.write(String.format("Variant: %s-%d-%s-%s, TranscriptAnnotations: %s, Jannovar Gene Symbol: %s, Jannovar Effect: %s\n", chrom, pos, ref, alt, x, jannovarGeneSymbol, jannovarFirstEffect));

                    int isMatch = 0;
                    if (clinvarFirstEffect.equals(jannovarFirstEffect) || (clinvarFirstEffect.contains("intron_variant") && jannovarFirstEffect.contains("intron_variant"))) {
                        isMatch = 1;
                    }

                    bw.write(String.format("%s-%d-%s-%s\t%s\t%s\t%s\t%s\t%d\n", chrom, pos, ref, alt, clinvarGeneSymbol, clinvarFirstEffect, jannovarGeneSymbol, jannovarFirstEffect, isMatch));
                }
            }

            for (String effect : uniqueEffects) {
                effectWriter.write(effect);
                effectWriter.newLine();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void testMethodsStart() {
        Path path = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Map<String, Integer> countMultipleMConGene = new HashMap<>();
        try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(path))))) {
            for (String line; (line = bufferedReader.readLine()) != null; ) {
                if (line.startsWith("#")) {
                    continue;
                }
//                    System.out.println(line);
                // get tokens by splitting lines
                String[] token = line.split("\t");
                var info = token[7];
                String geneInfo = "";
                String mc = "";
                // split by ; to get array of infoFields
                var infoTokens = info.split(";");
                var infoMap = new HashMap<String, String>();
                for (String infoToken : infoTokens) {
                    // put each info field split by (=) key and value into map
                    var x = infoToken.split("=");
                    var key = x[0];
                    var value = x[1];
                    if (key != null && value != null && value.contains(",")) {
                        infoMap.put(key, value);
                    }
                }
                System.out.printf("%s %s%n", infoMap.get("GENEINFO"), infoMap.get("MC"));


            }
        } catch (IOException e) {

        }

    }

    @Test
    void geneCountList() {
        Path path = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Map<String, Integer> geneCounts = new HashMap<>();

        try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(path))))) {
            for (String line; (line = bufferedReader.readLine()) != null; ) {
                if (line.startsWith("#")) {
                    continue;
                }
                String[] tokens = line.split("\t");
                String info = tokens[7];
                String[] infoTokens = info.split(";");

                Map<String, String> infoMap = new HashMap<>();
                for (String infoToken : infoTokens) {
                    String[] keyValue = infoToken.split("=");
                    if (keyValue.length == 2) {
                        String key = keyValue[0];
                        String value = keyValue[1];
                        infoMap.put(key, value);
                    }
                }

                String mc = infoMap.get("MC");
                if (mc != null && mc.contains(",")) {
                    String geneInfo = infoMap.get("GENEINFO");
                    if (geneInfo != null) {
                        String[] genePairs = geneInfo.split("\\|");
                        for (String genePair : genePairs) {
                            String geneSymbol = genePair.split(":")[0];
                            geneCounts.put(geneSymbol, geneCounts.getOrDefault(geneSymbol, 0) + 1);
                        }
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        geneCounts.forEach((k, v) -> System.out.println(k + ": " + v));
    }

    @Test
    void geneCountListtxt() {
        Path path = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Map<String, Integer> geneCounts = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(path))));
             BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/geneCountsNoComma.txt"))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] tokens = line.split("\t"), infoTokens = tokens[7].split(";");
                Map<String, String> infoMap = new HashMap<>();
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }
                String mc = infoMap.get("MC"), geneInfo = infoMap.get("GENEINFO");
                if (mc != null && geneInfo != null) { // && mc.contains(",")
                    for (String gp : geneInfo.split("\\|")) {
                        String gs = gp.split(":")[0];
                        geneCounts.put(gs, geneCounts.getOrDefault(gs, 0) + 1);
                    }
                }
            }
            geneCounts.forEach((k, v) -> {
                try {
                    bw.write(k + ": " + v);
                    bw.newLine();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void singleAndMultiMCCountPerGeneTxt() {
//        Path path = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Path path = Path.of("/Users/carlo/Desktop/jannovar/jannovar-cli-0.38/jannovar_annotated_grch37.vcf");
        Map<String, List<Integer>> geneCounts = new HashMap<>();
//        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(path))));
        try (BufferedReader br = new BufferedReader(new InputStreamReader(Files.newInputStream(path)));
             BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/JannovarGeneCountsCombinedAndTotal.txt"))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith("#")) continue;
                String[] tokens = line.split("\t"), infoTokens = tokens[7].split(";");
                Map<String, String> infoMap = new HashMap<>();
                for (String t : infoTokens) {
                    String[] kv = t.split("=");
                    if (kv.length == 2) infoMap.put(kv[0], kv[1]);
                }
                String mc = infoMap.get("MC"), geneInfo = infoMap.get("GENEINFO");
                if (mc != null && geneInfo != null) {
                    for (String geneString : geneInfo.split("\\|")) {
                        String geneSymbol = geneString.split(":")[0];
                        geneCounts.putIfAbsent(geneSymbol, new ArrayList<>(List.of(0, 0, 0)));
                        if (mc.contains(",")) {
                            geneCounts.get(geneSymbol).set(1, geneCounts.get(geneSymbol).get(1) + 1);
                        } else {
                            geneCounts.get(geneSymbol).set(0, geneCounts.get(geneSymbol).get(0) + 1);
                        }
                        geneCounts.get(geneSymbol).set(2, geneCounts.get(geneSymbol).get(2) + 1);
                    }
                }
            }
            bw.write("gene_symbol\tsingle_count\tmulti_count\ttotal\n");
            geneCounts.forEach((k, v) -> {
                try {
                    bw.write(k + "\t" + v.get(0) + "\t" + v.get(1) + "\t" + v.get(2));
                    bw.newLine();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    @Test
    void testTSV() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Path outputPath = Paths.get("/Users/carlo/Desktop/GeneTableStart/all_variants.tsv");

        Map<String, Map<String, Integer>> geneEffectCounts = new HashMap<>();
        Set<String> allEffects = new HashSet<>();

        try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))))) {
            for (String line; (line = bufferedReader.readLine()) != null; ) {
                if (line.startsWith("#")) {
                    continue;
                }

                String[] tokens = line.split("\t");
                String info = tokens[7];
                String[] infoTokens = info.split(";");

                Map<String, String> infoMap = new HashMap<>();
                for (String infoToken : infoTokens) {
                    String[] keyValue = infoToken.split("=");
                    if (keyValue.length == 2) {
                        infoMap.put(keyValue[0], keyValue[1]);
                    }
                }

                String mc = infoMap.get("MC");
                System.out.println(mc);
                if (mc != null ) { // && mc.contains(",")
                    String geneInfo = infoMap.get("GENEINFO");
                    if (geneInfo != null) {
                        String[] effects = mc.split(",");
                        String geneSymbol = geneInfo.split(":")[0];

                        Map<String, Integer> effectCounts = geneEffectCounts.getOrDefault(geneSymbol, new HashMap<>());

                        for (String effect : effects) {
                            String[] effectParts = effect.split("\\|");
                            if (effectParts.length > 1) {
                                String variantEffect = effectParts[1]; // for header
                                allEffects.add(variantEffect);
                                effectCounts.put(variantEffect, effectCounts.getOrDefault(variantEffect, 0) + 1);
                            }

                            geneEffectCounts.put(geneSymbol, effectCounts);
                        }
                    }
                }
            }

            try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
                writer.write("Gene Symbol");
                for (String effect : allEffects) {
                    writer.write("\t" + effect);
                }
                writer.newLine();

                for (Map.Entry<String, Map<String, Integer>> entry : geneEffectCounts.entrySet()) {
                    String gene = entry.getKey();
                    Map<String, Integer> effects = entry.getValue();

                    writer.write(gene);
                    for (String effect : allEffects) {
                        writer.write("\t" + effects.getOrDefault(effect, 0));
                    }
                    writer.newLine();
                }
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    void testTSVclinsig() {
        Path inputPath = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Path outputPath = Paths.get("/Users/carlo/Desktop/GeneTableStart/with_clinsig.tsv");

        // Mapping of Gene -> Effect -> Clinical Significance -> Count
        Map<String, Map<String, Map<String, Integer>>> geneEffectClinSigCounts = new HashMap<>();
        Set<String> allEffects = new HashSet<>();
        Set<String> allClinSig = new HashSet<>(Arrays.asList("Benign", "Likely_benign", "Pathogenic", "Likely_pathogenic", "Uncertain_significance"));

        try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(inputPath))))) {
            for (String line; (line = bufferedReader.readLine()) != null; ) {
                if (line.startsWith("#")) {
                    continue;
                }

                String[] tokens = line.split("\t");
                String info = tokens[7];
                String[] infoTokens = info.split(";");

                Map<String, String> infoMap = new HashMap<>();
                for (String infoToken : infoTokens) {
                    String[] keyValue = infoToken.split("=");
                    if (keyValue.length == 2) {
                        infoMap.put(keyValue[0], keyValue[1]);
                    }
                }

                String mc = infoMap.get("MC");
                if (mc != null) { // && mc.contains(",")
                    String geneInfo = infoMap.get("GENEINFO");
                    String clnsig = infoMap.getOrDefault("CLNSIG", null);
                    if (geneInfo != null) {
                        String[] effects = mc.split(",");
                        String geneSymbol = geneInfo.split(":")[0];
                        Map<String, Map<String, Integer>> effectClinSigMap = geneEffectClinSigCounts.getOrDefault(geneSymbol, new HashMap<>());

                        for (String effect : effects) {
                            String[] effectParts = effect.split("\\|");
                            if (effectParts.length > 1) {
                                String variantEffect = effectParts[1]; // for header
                                allEffects.add(variantEffect);
                                Map<String, Integer> clinSigCount = effectClinSigMap.getOrDefault(variantEffect, new HashMap<>());
                                clinSigCount.put(clnsig, clinSigCount.getOrDefault(clnsig, 0) + 1);
                                effectClinSigMap.put(variantEffect, clinSigCount);
                            }
                        }
                        geneEffectClinSigCounts.put(geneSymbol, effectClinSigMap);
                    }
                }
            }

            try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
                writer.write("Gene Symbol");
                for (String effect : allEffects) {
                    for (String clinSig : allClinSig) {
                        writer.write("\t" + effect + " " + clinSig);
                    }
                }
                writer.newLine();

                for (Map.Entry<String, Map<String, Map<String, Integer>>> geneEntry : geneEffectClinSigCounts.entrySet()) {
                    String gene = geneEntry.getKey();
                    writer.write(gene);
                    for (String effect : allEffects) {
                        for (String clinSig : allClinSig) {
                            Map<String, Integer> effectEntry = geneEntry.getValue().get(effect);
                            Integer count = 0;
                            if (effectEntry != null && effectEntry.containsKey(clinSig)) {
                                count = effectEntry.get(clinSig);
                            }
                            writer.write("\t" + count);
                        }
                    }
                    writer.newLine();
                }

            } catch (IOException e) {
                e.printStackTrace();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }




}




// count number of commas
// or split into array and count number of strings in array to get the
// or do object
// however count the appearance of those things where like more than 1 per gene and see than what it is and more ...
