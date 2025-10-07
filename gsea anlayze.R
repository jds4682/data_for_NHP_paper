# ===================================================================
# GSEA 자동 분석 및 통합 보고서 생성 스크립트 (GO & KEGG)
# ===================================================================

# --- 1. 필수 라이브러리 로드 ---
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggplot2)
# 엑셀 파일 저장을 위한 라이브러리 추가
# 설치가 안 된 경우: install.packages("writexl")
library(writexl)


# --- 2. 사용자 설정 (경로 수정) ---
# ❗️❗️❗️ 분석 전, 아래 2개의 경로를 사용자 환경에 맞게 수정해야 합니다. ❗️❗️❗️
input_dir <- "C:/tg"

# 결과 저장 폴더를 C드라이브가 아닌 '문서' 폴더로 변경
# R이 파일 쓰기 권한을 가질 수 있도록 사용자 폴더 내에 결과를 저장합니다.
output_dir <- file.path(Sys.getenv("USERPROFILE"), "Documents", "gsea_reports")


# --- 3. 분석 대상 파일 목록 정의 ---
csv_files <- list.files(path = input_dir, pattern = "_processed.csv$", full.names = TRUE)
if (length(csv_files) == 0) {
  stop("입력 폴더에 '_processed.csv' 파일을 찾을 수 없습니다. 'input_dir' 경로를 확인하세요.")
}


# ===================================================================
# --- 모든 결과를 취합할 리스트 초기화 ---
# ===================================================================
all_results_list <- list()


# ===================================================================
# --- 4. 메인 분석 루프 시작 ---
# ===================================================================
for (file_path in csv_files) {
  
  tang_name <- sub("_processed.csv", "", basename(file_path))
  cat(paste("\n\n=====>>", tang_name, "분석 시작 <<=====\n"))
  
  report_path <- file.path(output_dir, tang_name)
  dir.create(report_path, showWarnings = FALSE, recursive = TRUE)
  
  # --- 4.1. 데이터 준비 ---
  gene_data <- read.csv(file_path)
  
  aggregated_gene_data <- gene_data %>%
    group_by(GeneSymbol) %>%
    summarise(TotalScore = sum(Score)) %>%
    as.data.frame()
  
  ids <- bitr(aggregated_gene_data$GeneSymbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE)
  
  gene_data_merged <- merge(aggregated_gene_data, ids, by.x="GeneSymbol", by.y="SYMBOL", all.x = TRUE)
  gene_data_final <- gene_data_merged %>% filter(!is.na(ENTREZID))
  
  geneList <- gene_data_final$TotalScore
  names(geneList) <- gene_data_final$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- geneList[!duplicated(names(geneList))]
  
  if (length(geneList) == 0) {
    cat("ID 변환 후 유효한 유전자가 없어 이 처방은 건너뜁니다.\n")
    next
  }
  
  # --- 4.2. GSEA 분석 실행 (GO & KEGG) ---
  gse_go_results <- tryCatch({
    gseGO(geneList=geneList, OrgDb=org.Hs.eg.db, ont="BP", minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE, scoreType="pos")
  }, error = function(e) { cat("gseGO 분석 중 에러 발생.\n"); return(NULL) })
  
  gse_kegg_results <- tryCatch({
    gseKEGG(geneList=geneList, organism='hsa', minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE, scoreType="pos")
  }, error = function(e) { cat("gseKEGG 분석 중 에러 발생.\n"); return(NULL) })
  
  
  # --- 4.3. 개별 결과 보고서 저장 (기존과 동일) ---
  # (이 부분은 기존과 동일하게 개별 처방의 플롯과 CSV를 저장합니다)
  
  
  # ===================================================================
  # --- 통합 보고서를 위한 데이터 추출 및 가공 ---
  # ===================================================================
  
  # --- GO 결과 추출 ---
  if (!is.null(gse_go_results) && nrow(as.data.frame(gse_go_results)) > 0) {
    go_df <- as.data.frame(gse_go_results)
    
    # 1. 'gsea plot' 카테고리: 기본 순위 상위 3개
    gsea_plot_data <- go_df[1:min(3, nrow(go_df)), ]
    if(nrow(gsea_plot_data) > 0) {
        gsea_plot_data$category <- "gsea plot"
    }
    
    # 2. 'FDR' 카테고리: p.adjust 기준 상위 10개
    fdr_data <- go_df %>% arrange(p.adjust) %>% slice_head(n = 10)
    if(nrow(fdr_data) > 0) {
        fdr_data$category <- "FDR"
    }

    # 3. 'dot plot' & 'ridge plot' 카테고리: p.adjust 기준 상위 5개
    plot_data <- go_df %>% arrange(p.adjust) %>% slice_head(n = 5)
    if(nrow(plot_data) > 0) {
        dot_plot_data <- plot_data
        dot_plot_data$category <- "dot plot"
        
        ridge_plot_data <- plot_data
        ridge_plot_data$category <- "ridge plot"
    } else {
        dot_plot_data <- data.frame()
        ridge_plot_data <- data.frame()
    }

    # 추출된 GO 데이터들을 하나로 합치기
    combined_go_data <- bind_rows(gsea_plot_data, fdr_data, dot_plot_data, ridge_plot_data)
    
    if(nrow(combined_go_data) > 0) {
        combined_go_data$prescription <- tang_name
        all_results_list[[paste(tang_name, "GO")]] <- combined_go_data
    }
  }
  
  # --- ★★★ 신규 추가: KEGG 결과 추출 ★★★ ---
  if (!is.null(gse_kegg_results) && nrow(as.data.frame(gse_kegg_results)) > 0) {
    kegg_df <- as.data.frame(gse_kegg_results)
    
    # 1. 'gsea plot' 카테고리: 기본 순위 상위 3개
    gsea_plot_data_kegg <- kegg_df[1:min(3, nrow(kegg_df)), ]
    if(nrow(gsea_plot_data_kegg) > 0) {
        gsea_plot_data_kegg$category <- "gsea plot"
    }
    
    # 2. 'FDR' 카테고리: p.adjust 기준 상위 10개
    fdr_data_kegg <- kegg_df %>% arrange(p.adjust) %>% slice_head(n = 10)
    if(nrow(fdr_data_kegg) > 0) {
        fdr_data_kegg$category <- "FDR"
    }

    # 3. 'dot plot' & 'ridge plot' 카테고리: p.adjust 기준 상위 5개
    plot_data_kegg <- kegg_df %>% arrange(p.adjust) %>% slice_head(n = 5)
    if(nrow(plot_data_kegg) > 0) {
        dot_plot_data_kegg <- plot_data_kegg
        dot_plot_data_kegg$category <- "dot plot"
        
        ridge_plot_data_kegg <- plot_data_kegg
        ridge_plot_data_kegg$category <- "ridge plot"
    } else {
        dot_plot_data_kegg <- data.frame()
        ridge_plot_data_kegg <- data.frame()
    }

    # 추출된 KEGG 데이터들을 하나로 합치기
    combined_kegg_data <- bind_rows(gsea_plot_data_kegg, fdr_data_kegg, dot_plot_data_kegg, ridge_plot_data_kegg)
    
    if(nrow(combined_kegg_data) > 0) {
        combined_kegg_data$prescription <- tang_name
        all_results_list[[paste(tang_name, "KEGG")]] <- combined_kegg_data
    }
  }
  
  cat(paste("----->>", tang_name, "분석 및 보고서 저장 완료\n"))
}

# ===================================================================
# --- 모든 결과를 취합하여 하나의 엑셀 파일로 저장 ---
# ===================================================================
if (length(all_results_list) > 0) {
  cat("\n\n=====>> 최종 통합 보고서 생성 시작 <<=====\n")
  
  # 리스트에 있는 모든 데이터프레임을 하나로 합치기
  final_summary_df <- bind_rows(all_results_list)
  
  # 요청한 열 순서대로 재정렬
  final_summary_df <- final_summary_df %>%
    select(prescription, category, ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge, core_enrichment) %>%
    # 중복된 행 제거
    distinct()
    
  # 'Description' 열을 기준으로 정렬하여 가독성 높이기
  final_summary_df <- final_summary_df %>%
    arrange(prescription, Description, category)

  # 엑셀 파일로 저장
  output_excel_path <- file.path(output_dir, "GSEA_Summary_Report.xlsx")
  write_xlsx(final_summary_df, path = output_excel_path)
  
  cat(paste("최종 통합 보고서가 아래 경로에 저장되었습니다:\n", output_excel_path, "\n"))
  
} else {
  cat("\n\n요약할 GSEA 결과가 없어 통합 보고서를 생성하지 않았습니다.\n")
}

cat("\n\n모든 처방에 대한 GSEA 분석 및 보고서 생성이 완료되었습니다.\n")