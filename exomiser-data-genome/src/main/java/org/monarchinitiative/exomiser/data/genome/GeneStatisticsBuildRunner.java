package org.monarchinitiative.exomiser.data.genome;

import de.charite.compbio.jannovar.annotation.VariantEffect;
import org.h2.mvstore.MVMap;
import org.h2.mvstore.MVStore;
import org.h2.mvstore.MVStoreTool;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarGeneStatistics;
import org.monarchinitiative.exomiser.data.genome.model.BuildInfo;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import java.nio.file.Path;
import java.util.Map;

public class GeneStatisticsBuildRunner {

    private static final Logger logger = LoggerFactory.getLogger(GeneStatisticsBuildRunner.class);

    private final Path buildPath;
    private final BuildInfo buildInfo;
    private final ClinVarGeneStatistics clinVarGeneStatistics;

    public GeneStatisticsBuildRunner(BuildInfo buildInfo, Path buildPath, ClinVarGeneStatistics clinVarGeneStatistics) {
        this.buildPath = buildPath;
        this.buildInfo = buildInfo;
        this.clinVarGeneStatistics = clinVarGeneStatistics;
    }

    public void run() {
        String fileName = buildPath.resolve(buildInfo.getBuildString() + "_genestats.mv.db").toString();
        MVStore mvStore = new MVStore.Builder()
                .fileName(fileName)
                .compress()
                .open();

        mvStore.setVersionsToKeep(0);

        MVMap<String, Map<VariantEffect, Map<ClinVarData.ClinSig, Integer>>> geneStatsMap = mvStore.openMap("GeneStatistics");

        geneStatsMap.putAll(clinVarGeneStatistics.getMap());

        logger.info("Written {} gene statistics to store", geneStatsMap.size());

        mvStore.commit();
        mvStore.compactMoveChunks();

        logger.info("Compacting store...");
        MVStoreTool.compact(fileName, true);
    }
}
