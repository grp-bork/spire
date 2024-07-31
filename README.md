# SPIRE workflow
<table>
  <tr width="100%">
    <td width="150px">
      <a href="https://www.bork.embl.de/"><img src="https://www.bork.embl.de/assets/img/normal_version.png" alt="Bork Group Logo" width="150px" height="auto"></a>
    </td>
    <td width="425px" align="center">
      <b>Developed by the <a href="https://www.bork.embl.de/">Bork Group</a></b><br>
      Raise an <a href="https://github.com/grp-bork/spire/issues">issue</a> or <a href="mailto:N4M@embl.de">contact us</a><br><br>
      See our <a href="https://www.bork.embl.de/services.html">other Software & Services</a>
    </td>
    <td width="500px">
      Contributors:<br>
      <ul>
        <li>
  <a href="https://github.com/fullama/">Anthony Fullam</a> <a href="https://orcid.org/0000-0002-0884-8124"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
</li>
        <li>
  <a href="https://github.com/defleury/">Thomas Sebastian B. Schmidt</a> <a href="https://orcid.org/0000-0001-8587-4177"><img src="https://orcid.org/assets/vectors/orcid.logo.icon.svg" alt="ORCID icon" width="20px" height="20px"></a><br>
</li>
      </ul>
    </td>
  </tr>
</table>

</table>

---
#### Description

#### Citation

Cite the [SPIRE publication](https://doi.org/10.1093/nar/gkad943) when using our workflow:
```
Schmidt TSB, Fullam A, Ferretti P, et al. SPIRE: a Searchable, Planetary-scale mIcrobiome REsource. Nucleic Acids Res. 2024;52(D1):D777-D783. doi:10.1093/nar/gkad943
```

---
# Overview

**Preprocessing:**
1. Trimming (ngless)
2. Length filtering (ngless)
3. Human DNA decontamination (GRCh38)

**MAGs:**
1. Assembly (megahit)
2. Gene Calling (prodigal)
3. Remove Small Contigs (seqtk)
4. Index (bwa)
5. Alignment (bwa, samtools)
6. Calculate Depths (jgi_summarize_bam_contig_depths)
7. Binning (metabat2)
8. Per-Bin Gene Calling (seqtk)
9. Assembly Stats (assembly-stats)
10. Assembly Mash Sketching (mash)
11. Bin Mash Sketching (mash)

**Annotation:**
1. rRNA Detection (barrnap)
2. ARG Annotation (abricate, rgi)
3. Virulence Factor Annotation (abricate)
4. sORFs Detection (macrel)
5. Genome Quality Assessment (checkm2, gunc)
6. Functional Annotation (eggnog-mapper)
7. Taxonomic Classification (gtdbtk)
