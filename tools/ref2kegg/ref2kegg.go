package main

import (
	//"compress/gzip"
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	//"reflect"
	"regexp"
	//"time"
	//"runtime"
	"strconv"
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
	r, _ := regexp.Compile("^>(\\S+) ")
	childr, _ := regexp.Compile("^(\\S+)")
	for {
		linestr, err := reader.ReadString('\n')
		line := []byte(linestr)
		if err == io.EOF {
			break
		}
		isHeader, err := regexp.Match(">", line[:1])
		if err != nil {
			panic(err)
		}
		if isHeader {
			parent := r.FindAllStringSubmatch(string(line), -1)[0][1]
			//fmt.Println(string(line))
			//fmt.Println(parent)
			for i, rune := range bytes.SplitN(line, []byte{1}, -1) {
				if i != 0 {
					children := childr.FindAllString(string(rune), -1)

					if len(children) > 0 {
						child := children[0]
						//fmt.Printf("counter %d: %s\n", i, child)
						m[child] = parent
						//fmt.Printf("Acc to WP - acc:%s\trow:%s\n", child, parent)
					} else {
						log.Println("Problematic:")
						//fmt.Println(string(line))
						//fmt.Printf("%s\n", string(rune))
					}
				}
			}
		}
	}
	//for key, _ := range m {
	//fmt.Println(reflect.TypeOf(key))
	//}
	//if "one" == "one" {

	//fmt.Println("test")
	//}
	//fmt.Println(reflect.TypeOf("NP_065910"))
	return m
}

func accKEGG(filepath string) map[string]string {
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
	return m
}

func KEGGko(filepath string) map[string]int {
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

	//start := time.Now()

	flag.Parse()
	nrpath := flag.Arg(0)
	linker_keggGene2acc := flag.Arg(1)
	linker_koKeggGene := flag.Arg(2)

	//accNRACC_map := accNRACC(nrpath)
	//for stuff, element := range accNRACC_map {
	//fmt.Printf("%s\t%s\n", element, stuff)
	//}

	accNRACC_map := accNRACC(nrpath)
	accKEGG_map := accKEGG(linker_keggGene2acc)
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
	//elapsed := time.Since(start)
	//log.Printf("Binomial took %s", elapsed)
}
