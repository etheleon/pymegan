package main

import (
	//"compress/gzip"
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"regexp"
	"time"
	//"runtime"
	//"strconv"
	"strings"
	//"sync"
)

//example
//>WP_003131952.1 30S ribosomal protein S18 [Lactococcus lactis]^ANP_268346.1 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis Il1403]^AQ9CDN0.1 RecName: Full=30S ribosomal protein S18^AQ02VU1.1 RecName: Full=30S ribosomal protein S18^AA2RNZ2.1 RecName: Full=30S ribosomal protein S18^AAAK06287.1 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis Il1403]^AABJ73931.1 SSU ribosomal protein S18P [Lactococcus lactis subsp. cremoris SK11]^ACAL99037.1 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris MG1363]^AADA65983.1 SSU ribosomal protein S18P [Lactococcus lactis subsp. lactis KF147]^AADJ61439.1 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris NZ9000]^AADZ64834.1 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis CV56]^AEHE92602.1 hypothetical protein LLCRE1631_01913 [Lactococcus lactis subsp. lactis CNCM I-1631]

func accNRACC(nrpath string) map[string]string {
	m := make(map[string]string)
	//open file
	nr, err := os.Open(nrpath)
	if err != nil {
		log.Fatal(err)
	}
	defer nr.Close()

	reader := bufio.NewReader(nr)
	r, _ := regexp.Compile("^>(\\W+) ")
	childr, _ := regexp.Compile("^(\\W+) ")
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		isHeader, err := regexp.Match(">", line[:1])
		if err != nil {
			panic(err)
		}
		if isHeader {
			parent := r.FindAllString(string(line), -1)[0]
			children := strings.Split(line, "\x01")
			for _, element := range children {
				childACC := childr.FindAllString(string(element), -1)
				m[childACC] = parent
			}
		}
	}
	return m
}

func main() {

	start := time.Now()

	flag.Parse()
	nrpath := flag.Arg(0)
	accNRACC(nrpath)
	//linker_keggGene2gi := flag.Arg(1)
	//linker_koKeggGene := flag.Arg(2)

	//giWPGI_map := giWPGI(nrpath)
	//giKEGG_map := giKEGG(linker_keggGene2gi)
	//KEGGko_map := KEGGko(linker_koKeggGene)

	//Link 3 dicts (KEGGko::ko <- KEGG geneID) <-> (kegggi::KEGG geneID <- gi) <-> (giwp::gi <- WP_refseq's gi)

	//    keggGI_map := make(map[string]int)
	//    for gi, geneID := range giKEGG_map {
	//        //fmt.Printf("%s\t%d\n", geneID, giWPGI_map[gi])
	//        if WPGI, ok := giWPGI_map[gi]; ok {
	//            keggGI_map[geneID] = WPGI
	//        }
	//    }
	//
	//    wpgiKO_map := make(map[int]int)
	//
	//    for geneID, ko := range KEGGko_map {
	//        //fmt.Printf("%s\t%d\n", geneID, keggGI_map[geneID])
	//        if gi, ok := keggGI_map[geneID]; ok {
	//            wpgiKO_map[gi] = ko
	//        }
	//    }
	//    for wpgi, ko := range wpgiKO_map {
	//        fmt.Printf("%d\t%d\n", wpgi, ko)
	//    }
	elapsed := time.Since(start)
	log.Printf("Binomial took %s", elapsed)

}
