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
	"strconv"
	"strings"
)

func giWP(nrpath string) map[int]string {
	//m(gi) = WP (non-redundant WP)
	m := make(map[int]string)
	nr, err := os.Open(nrpath)
	if err != nil {
		log.Fatal(err)
	}
	defer nr.Close()

	reader := bufio.NewReader(nr)
	r, _ := regexp.Compile("gi\\|\\d+\\|ref\\|[^|]+")
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		//header := string(line[:1])
		isHeader, err := regexp.Match(">", line[:1])
		if err != nil {
			panic(err)
		}
		if isHeader {
			matches := r.FindAllString(string(line), -1)
			switch {
			case len(matches) > 1:
				nrEntry, children := matches[0], matches[1:]

				parent := strings.Split(nrEntry, "|")
				parentRefSeq := parent[3]
				giInt, _ := strconv.Atoi(parent[1])
				//fmt.Println(giInt)
				//fmt.Println(parentRefSeq)
				m[giInt] = parentRefSeq
				for _, element := range children {
					results := strings.Split(element, "|")
					gi, _ := strconv.Atoi(results[1])
					//childRefSeq := results[3]
					m[gi] = parentRefSeq
					//fmt.Printf("%s\t%s\t%s\n", parentRefSeq, gi, childRefSeq)
				}
			case len(matches) == 1:
				parent := strings.Split(matches[0], "|")
				parentRefSeq := parent[3]
				giInt, _ := strconv.Atoi(parent[1])
				m[giInt] = parentRefSeq
				//case len(matches) == 0:
				//fmt.Println(matches)
			}
		}
	}
	return m
}

func giKEGG(filepath string) map[int]string {
	m := make(map[int]string)
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatal(err)
	}
	giRegex, _ := regexp.Compile("ncbi-gi:(\\d+)")

	defer file.Close()
	reader := bufio.NewReader(file)
	for {
		line, _, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		//log.Println(scanner.Text())
		row := strings.Split(string(line), "\t")
		matches := giRegex.FindStringSubmatch(row[1])
		//fmt.Println(matches[1])
		gi, _ := strconv.Atoi(matches[1])
		m[gi] = row[0]
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
		//log.Println(scanner.Text())
		row := strings.Split(string(line), "\t")

		matches := koRegex.FindStringSubmatch(row[1])
		//fmt.Println(matches[1])
		ko, _ := strconv.Atoi(matches[1])
		m[matches[0]] = ko
	}
	return m
}

func main() {
	flag.Parse()
	nrpath := flag.Arg(0)
	linker_keggGene2gi := flag.Arg(1)
	linker_koKeggGene := flag.Arg(2)

	//giWP(nrpath)
	giWP_map := giWP(nrpath)
	giKEGG_map := giKEGG(linker_keggGene2gi)
	KEGGko_map := KEGGko(linker_koKeggGene)

	//Link 3 dicts (ko <- KEGG geneID) <-> (kegggi::KEGG geneID <- gi) <-> (giwp::gi <- WP_refseq)
	var keggWP_map map[string]string

	for gi, geneID := range giKEGG_map {
		WP := giWP_map[gi]
		keggWP_map[geneID] = WP
	}

	var wpKO_map map[string]int

	for geneID, ko := range KEGGko_map {
		WP := keggWP_map[geneID]
		wpKO_map[WP] = ko
	}

	for wp, ko := range wpKO_map {
		fmt.Printf("%s\t%s\n", wp, ko)
	}
}
