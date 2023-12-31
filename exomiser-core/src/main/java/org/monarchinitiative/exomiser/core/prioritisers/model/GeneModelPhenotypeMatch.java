/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2021 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.monarchinitiative.exomiser.core.prioritisers.model;

import com.fasterxml.jackson.annotation.JsonIgnore;
import org.monarchinitiative.exomiser.core.phenotype.ModelPhenotypeMatch;
import org.monarchinitiative.exomiser.core.phenotype.Organism;
import org.monarchinitiative.exomiser.core.phenotype.PhenotypeMatch;

import java.util.List;
import java.util.Objects;

/**
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class GeneModelPhenotypeMatch implements GeneModel {

    private final double score;
    private final GeneModel model;
    private final List<PhenotypeMatch> bestModelPhenotypeMatches;

    public GeneModelPhenotypeMatch(double score, GeneModel model, List<PhenotypeMatch> bestModelPhenotypeMatches) {
        this.score = score;
        this.model = model;
        this.bestModelPhenotypeMatches = List.copyOf(bestModelPhenotypeMatches);
    }

    public GeneModelPhenotypeMatch(ModelPhenotypeMatch<GeneModel> match) {
        this(match.getScore(), match.getModel(), match.getBestPhenotypeMatches());
    }

    public double getScore() {
        return score;
    }

    public GeneModel getModel() {
        return model;
    }

    public List<PhenotypeMatch> getBestModelPhenotypeMatches() {
        return bestModelPhenotypeMatches;
    }

    @JsonIgnore
    @Override
    public String getId() {
        return model.getId();
    }

    @JsonIgnore
    @Override
    public List<String> getPhenotypeIds() {
        return model.getPhenotypeIds();
    }

    @JsonIgnore
    @Override
    public Organism getOrganism() {
        return model.getOrganism();
    }

    @JsonIgnore
    @Override
    public int getEntrezGeneId() {
        return model.getEntrezGeneId();
    }

    @JsonIgnore
    @Override
    public String getHumanGeneSymbol() {
        return model.getHumanGeneSymbol();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GeneModelPhenotypeMatch)) return false;
        GeneModelPhenotypeMatch that = (GeneModelPhenotypeMatch) o;
        return Double.compare(that.score, score) == 0 &&
                Objects.equals(model, that.model) &&
                Objects.equals(bestModelPhenotypeMatches, that.bestModelPhenotypeMatches);
    }

    @Override
    public int hashCode() {
        return Objects.hash(model, bestModelPhenotypeMatches, score);
    }


    @Override
    public String toString() {
        return "GeneModelPhenotypeMatch{" +
                "score=" + score +
                ", model=" + model +
                ", bestModelPhenotypeMatches=" + bestModelPhenotypeMatches +
                '}';
    }
}
