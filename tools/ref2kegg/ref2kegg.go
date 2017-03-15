package main

import (
	//"compress/gzip"
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	//"reflect"
	"regexp"
	"time"
	//"runtime"
	"strconv"
	"strings"
	//"sync"
)

//example
//>WP_003131952.1 30S ribosomal protein S18 [Lactococcus lactis]^ANP_268346.1 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis Il1403]^AQ9CDN0.1 RecName: Full=30S ribosomal protein S18^AQ02VU1.1 RecName: Full=30S ribosomal protein S18^AA2RNZ2.1 RecName: Full=30S ribosomal protein S18^AAAK06287.1 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis Il1403]^AABJ73931.1 SSU ribosomal protein S18P [Lactococcus lactis subsp. cremoris SK11]^ACAL99037.1 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris MG1363]^AADA65983.1 SSU ribosomal protein S18P [Lactococcus lactis subsp. lactis KF147]^AADJ61439.1 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris NZ9000]^AADZ64834.1 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis CV56]^AEHE92602.1 hypothetical protein LLCRE1631_01913 [Lactococcus lactis subsp. lactis CNCM I-1631]

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func accNRACC(nrpath string, legit map[string]string) map[string]string {
	m := make(map[string]string)
	r, _ := regexp.Compile("^>([^\\s\\.]+)")
	childr, _ := regexp.Compile("^[^\\s\\.]+")

	//slurp file
	log.Printf("Slurping file: %s\n", nrpath)
	dat, err := ioutil.ReadFile(nrpath)
	check(err)
	log.Printf("Finished reading in: %s\n", nrpath)

	for _, line := range bytes.Split(dat, []byte{'\n'}) {
		//fmt.Println(string(line))
		isHeader, err := regexp.Match(">", line[:1])
		check(err)
		if isHeader {
			parent := r.FindAllStringSubmatch(string(line), -1)[0][1]
			//fmt.Printf("%shehhlo\n", parent)
			for i, rune := range bytes.SplitN(line, []byte{1}, -1) {
				if i != 0 {
					children := childr.FindAllString(string(rune), -1)
					if len(children) > 0 {
						child := children[0]
						//log.Printf("Acc to WP - acc:%stest\trow:%s\n", child, parent)
						if _, ok := legit[child]; ok {
							//log.Printf("This acc:%s has a KEGG record", child)
							m[child] = parent
						}
						//fmt.Printf("Acc to WP - acc:%stest\trow:%s\n", child, parent)
					} else {
						log.Println("Problematic:")
					}
				}
			}
		}
	}
	return m
}

func accKEGG(filepath string) map[string]string {
	log.Println("Storing ncbi geneid to kegg geneid...")
	m := make(map[string]string)
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatal(err)
	}
	accRegex, _ := regexp.Compile("ncbi-proteinid:(\\S+)")

	defer file.Close()
	reader := bufio.NewReader(file)
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		row := strings.Split(string(line), "\t")
		matches := accRegex.FindStringSubmatch(row[1])
		acc := matches[1]
		m[acc] = row[0]
		//fmt.Printf("Acc to KEGG - acc:%s\tKEGGid:%s\n", acc, row[0])
	}
	log.Println("Finished storing kegg geneid to ncbi geneid")
	return m
}

func KEGGko(filepath string) map[string]int {
	log.Println("Storing KO to kegg GeneID...")
	m := make(map[string]int)
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatal(err)
	}
	koRegex, _ := regexp.Compile("ko:K(\\d+)")

	defer file.Close()
	reader := bufio.NewReader(file)
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		row := strings.Split(string(line), "\t")

		matches := koRegex.FindStringSubmatch(row[1])
		ko, _ := strconv.Atoi(matches[1])
		m[row[0]] = ko
		//fmt.Printf("KEGGid to KO - KEGGid:%s\tko:%d\n", row[0], ko)
	}
	return m
}

func main() {

	start := time.Now()

	flag.Parse()
	nrpath := flag.Arg(0)
	linker_keggGene2acc := flag.Arg(1)
	linker_koKeggGene := flag.Arg(2)

	//accNRACC_map := accNRACC(nrpath)
	//for stuff, element := range accNRACC_map {
	//fmt.Printf("%s\t%s\n", element, stuff)
	//}

	accKEGG_map := accKEGG(linker_keggGene2acc)
	accNRACC_map := accNRACC(nrpath, accKEGG_map)
	KEGGko_map := KEGGko(linker_koKeggGene)

	//Link 3 dicts (KEGGko::ko <- KEGG geneID) <-> (kegggi::KEGG geneID <- gi) <-> (giwp::gi <- WP_refseq's gi)
	keggAcc_map := make(map[string]string)
	for acc, geneID := range accKEGG_map {
		if nracc, ok := accNRACC_map[acc]; ok {
			keggAcc_map[geneID] = nracc
		}
	}

	nraccKO_map := make(map[string]int)

	for geneID, ko := range KEGGko_map {
		//fmt.Printf("%s\t%d\n", geneID, keggGI_map[geneID])
		if acc, ok := keggAcc_map[geneID]; ok {
			nraccKO_map[acc] = ko
		}
	}
	for nracc, ko := range nraccKO_map {
		fmt.Printf("%s\t%d\n", nracc, ko)
	}
	elapsed := time.Since(start)
	log.Printf("Binomial took %s", elapsed)
}
