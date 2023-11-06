package org.monarchinitiative.exomiser.data.genome.model.resource;


import org.junit.jupiter.api.Test;
import org.monarchinitiative.exomiser.core.genome.GenomeAssembly;
import org.monarchinitiative.exomiser.core.genome.JannovarVariantAnnotator;
import org.monarchinitiative.exomiser.core.genome.VariantAnnotator;
import org.monarchinitiative.exomiser.core.model.ChromosomalRegionIndex;
import org.monarchinitiative.exomiser.core.model.VariantAnnotation;
import org.monarchinitiative.svart.CoordinateSystem;
import org.monarchinitiative.svart.GenomicVariant;
import org.monarchinitiative.svart.Strand;
import org.monarchinitiative.svart.util.VariantTrimmer;


import java.awt.image.DataBufferDouble;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.zip.GZIPInputStream;



class ClinVarGeneStatsParserTest {
//    private final JannovarVariantAnnotator instance = new JannovarVariantAnnotator(TestFactory.getDefaultGenomeAssembly(), TestFactory
//            .buildDefaultJannovarData(), ChromosomalRegionIndex.empty());

    private List<VariantAnnotation> annotate(VariantAnnotator instance, String contig, int start, String ref, String alt) {
        GenomicVariant variant = variant(contig, start, ref, alt);
        return instance.annotate(variant);
    }

    private GenomicVariant variant(String contig, int start, String ref, String alt) {
        VariantTrimmer.VariantPosition variantPosition = VariantTrimmer.leftShiftingTrimmer(VariantTrimmer.retainingCommonBase()).trim(Strand.POSITIVE, start, ref, alt);
        return GenomicVariant.of(GenomeAssembly.HG19.getContigByName(contig), Strand.POSITIVE, CoordinateSystem.ONE_BASED, variantPosition.start(), variantPosition.ref(), variantPosition.alt());
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
        Path path = Path.of("/Users/carlo/Downloads/clinvar.vcf.gz");
        Map<String, List<Integer>> geneCounts = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(path))));
             BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/carlo/Desktop/GeneTableStart/geneCountsCombinedAndTotal.tsv"))) {
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