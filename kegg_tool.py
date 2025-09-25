#!/usr/bin/env python3
"""
KEGG数据库访问工具
用于查询pathway、ko、genome、reaction、enzyme、genes等信息
(版本：最终修复版 - 解决初始化冲突)
"""

from crewai.tools import BaseTool
import requests
import json
from typing import Dict, List, Optional, Any, Union
from pydantic import BaseModel, Field, PrivateAttr
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import time
import argparse

# ---- KEGG flat-text 轻量解析器 (无变动) ----

def _parse_kegg_flat(text: str) -> Dict[str, Union[str, List[str]]]:
    """将 KEGG get 返回的平文本解析为一个 dict（保留常用字段）"""
    out: Dict[str, Union[str, List[str]]] = {}
    cur_key = None
    for raw in text.splitlines():
        if not raw.strip(): continue
        key = raw[:12].strip()
        val = raw[12:].strip()
        if key:
            cur_key = key
            if key not in out: out[key] = []
        if cur_key is None: continue
        assert isinstance(out[cur_key], list)
        out[cur_key].append(val)
    for k, v in list(out.items()):
        if isinstance(v, list):
            if len(v) == 1: out[k] = v[0]
            else: out[k] = v
    return out

def _parse_kegg_flat_batch(text: str) -> Dict[str, Dict[str, Any]]:
    """将 KEGG get batch 返回的平文本解析为多个 dict 的集合"""
    results = {}
    entries_text = text.strip().split('\n///\n')
    for entry_text in entries_text:
        if not entry_text.strip(): continue
        parsed_entry = _parse_kegg_flat(entry_text)
        entry_id_val = parsed_entry.get("ENTRY")
        if isinstance(entry_id_val, str):
            entry_id = entry_id_val.split()[0]
            results[entry_id] = parsed_entry
    return results

def _split_tsv(text: str) -> List[Dict[str, str]]:
    """将两列的TSV文本分割成对象列表"""
    items: List[Dict[str, str]] = []
    for line in text.strip().split('\n'):
        if not line: continue
        parts = line.split('\t')
        if len(parts) == 2: items.append({"id": parts[0], "value": parts[1]})
        elif len(parts) == 1: items.append({"id": parts[0], "value": ""})
    return items

# ---------------- 参数模型 (无变动) ----------------
class KeggToolInput(BaseModel):
    operation: Optional[str] = Field(None, description="要执行的操作，例如 'find_compound_by_name', 'resolve_compound_workflow'")
    database: Optional[str] = Field(None, description="KEGG 数据库名, 如 'pathway', 'compound', 'ko'")
    organism: Optional[str] = Field(None, description="物种代码, 如 'hsa' (人), 'eco' (大肠杆菌)")
    keywords: Optional[str] = Field(None, description="用于 find 操作的搜索关键词")
    entry_id: Optional[str] = Field(None, description="KEGG 条目ID, 如 'C00031', 'ko:K00844'")
    format_type: Optional[str] = Field("txt", description="返回格式 (txt, aaseq, ntseq等)")
    target_db: Optional[str] = Field(None, description="用于 link/convert 操作的目标数据库")
    source_db_entries: Optional[str] = Field(None, description="用于 link 操作的源数据库条目")
    source_ids: Optional[str] = Field(None, description="用于 convert 操作的源ID")
    compound_id: Optional[str] = Field(None, description="化合物ID, 如 'C00031'")
    pathway_id: Optional[str] = Field(None, description="通路ID, 如 'map00010'")
    name: Optional[str] = Field(None, description="化合物或基因等的名称")
    reaction_id: Optional[str] = Field(None, description="反应ID, 如 'R00225'")
    name_or_id: Optional[str] = Field(None, description="在 workflow 中使用的化合物名称或ID")
    max_items: Optional[int] = Field(20, description="在 workflow 中限制返回结果的数量")

# ---------------- 工具类 (结构性修复) ----------------
class KeggTool(BaseTool):
    name: str = "kegg_tool"
    description: str = "KEGG 查询工具，可以查找化合物、通路、反应等信息，并支持一键式工作流'resolve_compound_workflow'来获取化合物的完整信息。"
    args_schema = KeggToolInput

    # [结构性修复] 移除自定义 __init__，将配置改为类属性，以兼容 Pydantic/CrewAI 的初始化
    base_url: str = "https://rest.kegg.jp"
    timeout: float = 15.0
    max_retries: int = 3
    
    # 使用 PrivateAttr 来存储不属于模型公开接口的状态
    _session: Optional[requests.Session] = PrivateAttr(default=None)

    def _get_session(self) -> requests.Session:
        """延迟初始化并返回 requests.Session 对象"""
        if self._session is None:
            sess = requests.Session()
            adapter = HTTPAdapter(max_retries=Retry(total=self.max_retries, backoff_factor=0.3,
                                                   status_forcelist=[429, 500, 502, 503, 504]))
            sess.mount("http://", adapter)
            sess.mount("https://", adapter)
            sess.headers.update({"User-Agent": "Biocrew-KeggTool/1.0"})
            self._session = sess
        return self._session

    def _run(
        self,
        operation: Optional[str] = None, database: Optional[str] = None, organism: Optional[str] = None,
        keywords: Optional[str] = None, entry_id: Optional[str] = None, format_type: str = "txt",
        target_db: Optional[str] = None, source_db_entries: Optional[str] = None, source_ids: Optional[str] = None,
        compound_id: Optional[str] = None, pathway_id: Optional[str] = None, name: Optional[str] = None,
        reaction_id: Optional[str] = None, name_or_id: Optional[str] = None, max_items: int = 20,
    ) -> str:
        """CrewAI BaseTool entrypoint (dispatcher)."""
        op = (operation or "").strip().lower()
        if not op and (name_or_id or name): op = "resolve_compound_workflow"

        result_dict = {}
        try:
            if op == "find_compound_by_name":
                result_dict = self.find_compound_by_name(name or keywords or "")
            elif op in {"workflow", "resolve_compound_workflow"}:
                search_term = name_or_id or name or ""
                result_dict = self.resolve_compound_workflow(name_or_id=search_term, organism=organism, max_items=max_items)
            # ... (可以根据需要添加其他 op 的 elif 分支)
            else:
                result_dict = {"status": "error", "message": f"Unsupported or missing operation: {op or '[none]'}"}
        except Exception as e:
            result_dict = {"status": "error", "message": f"dispatcher error in '{op}': {type(e).__name__}: {e}"}
        
        return json.dumps(result_dict, ensure_ascii=False, indent=2)

    def run(self, **kwargs) -> str:
        """同步调用入口，兼容某些 BaseTool 版本要求存在 run 方法"""
        return self._run(**kwargs)

    async def _arun(self, **kwargs) -> str:
        """异步调用入口，兼容某些 BaseTool 版本要求存在 _arun 方法"""
        return self._run(**kwargs)

    def find_compound_by_name(self, name_or_keywords: str) -> Dict:
        try:
            if not name_or_keywords: return {"status": "error", "message": "缺少关键词"}
            kw = name_or_keywords.replace(' ', '+')
            url = f"{self.base_url}/find/compound/{kw}"
            r = self._get_session().get(url, timeout=self.timeout)
            r.raise_for_status()
            return {"status": "success", "data": _split_tsv(r.text)}
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def pathways_of_compound(self, compound_id: str) -> Dict:
        try:
            cid = compound_id if compound_id.lower().startswith('c') else compound_id
            url = f"{self.base_url}/link/pathway/{cid}"
            r = self._get_session().get(url, timeout=self.timeout); r.raise_for_status()
            items: List[Dict[str, str]] = []
            for line in r.text.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    left, right = parts[0], parts[1]
                    pid = left if left.startswith('path:') else (right if right.startswith('path:') else None)
                    if pid:
                        items.append({"id": pid, "value": compound_id})
                elif len(parts) == 1 and parts[0].startswith('path:'):
                    items.append({"id": parts[0], "value": compound_id})
            return {"status": "success", "data": items}
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def reactions_of_compound(self, compound_id: str) -> Dict:
        try:
            cid = compound_id if compound_id.lower().startswith('c') else compound_id
            url = f"{self.base_url}/link/reaction/{cid}"
            r = self._get_session().get(url, timeout=self.timeout); r.raise_for_status()
            items: List[Dict[str, str]] = []
            for line in r.text.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    left, right = parts[0], parts[1]
                    rid = left if left.startswith('rn:') else (right if right.startswith('rn:') else None)
                    if rid:
                        items.append({"id": rid, "value": compound_id})
                elif len(parts) == 1 and parts[0].startswith('rn:'):
                    items.append({"id": parts[0], "value": compound_id})
            return {"status": "success", "data": items}
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def enzymes_of_compound(self, compound_id: str) -> Dict:
        try:
            cid = compound_id if compound_id.lower().startswith('c') else compound_id
            url = f"{self.base_url}/link/enzyme/{cid}"
            r = self._get_session().get(url, timeout=self.timeout); r.raise_for_status()
            items: List[Dict[str, str]] = []
            for line in r.text.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    left, right = parts[0], parts[1]
                    eid = left if left.startswith('ec:') else (right if right.startswith('ec:') else None)
                    if eid:
                        items.append({"id": eid, "value": compound_id})
                elif len(parts) == 1 and parts[0].startswith('ec:'):
                    items.append({"id": parts[0], "value": compound_id})
            return {"status": "success", "data": items}
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def kos_of_reaction(self, reaction_id: str) -> Dict:
        try:
            rid = reaction_id if reaction_id.lower().startswith('rn:') else f"rn:{reaction_id}"
            url = f"{self.base_url}/link/ko/{rid}"
            r = self._get_session().get(url, timeout=self.timeout); r.raise_for_status()
            items: List[Dict[str, str]] = []
            for line in r.text.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) == 2:
                    left, right = parts[0], parts[1]
                    kid = left if left.startswith('ko:') else (right if right.startswith('ko:') else None)
                    if kid:
                        items.append({"id": kid, "value": rid})
                elif len(parts) == 1 and parts[0].startswith('ko:'):
                    items.append({"id": parts[0], "value": rid})
            return {"status": "success", "data": items}
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def genes_of_ko(self, ko_id: str, organism: Optional[str] = None) -> Dict:
        try:
            kid = ko_id if ko_id.lower().startswith('ko:') else f"ko:{ko_id}"
            url = f"{self.base_url}/link/genes/{kid}"
            r = self._get_session().get(url, timeout=self.timeout); r.raise_for_status()
            items = _split_tsv(r.text)
            if organism:
                prefix = f"{organism.lower()}:"
                items = [it for it in items if it['id'].lower().startswith(prefix) or it['value'].lower().startswith(prefix)]
            return {"status": "success", "data": items}
        except Exception as e:
            return {"status": "error", "message": str(e)}

    def resolve_compound_workflow(self, name_or_id: str, organism: Optional[str] = None, max_items: int = 20, return_envelope: bool = False) -> Dict[str, Any]:
        # 1) 解析/解析出 compound_id
        compound_id = None
        resolved_label = name_or_id
        if name_or_id.lower().startswith('c') and name_or_id[1:].isdigit():
            compound_id = name_or_id.upper()
        else:
            found = self.find_compound_by_name(name_or_id)
            if found.get('status') == 'success' and found.get('data'):
                # 优先精确匹配（按分号拆分同义词表）
                q = name_or_id.strip().lower()
                candidates = found['data']
                exact = None
                for it in candidates:
                    names = [s.strip().lower() for s in it['value'].split(';')]
                    if q in names:
                        exact = it
                        break
                chosen = exact or candidates[0]
                compound_id = chosen['id'].split(':')[-1]
                resolved_label = chosen['value'].split(';')[0]
            else:
                return {"status": "error", "message": f"未找到化合物 {name_or_id}"}

        # 2) 通路/反应/酶
        pw = self.pathways_of_compound(compound_id)
        pathways = [it.get('id') for it in pw.get('data', [])][:max_items]

        rx = self.reactions_of_compound(compound_id)
        reaction_ids = [it.get('id') for it in rx.get('data', [])][:max_items]

        # 批量取反应详情
        reactions: List[Dict[str, Any]] = []
        if reaction_ids:
            try:
                batch_ids = "+".join(reaction_ids)
                url = f"{self.base_url}/get/{batch_ids}"
                r = self._get_session().get(url, timeout=self.timeout); r.raise_for_status()
                batch = _parse_kegg_flat_batch(r.text)
                for rid_full in reaction_ids:
                    rid_short = rid_full.split(':')[-1]
                    detail = batch.get(rid_short, {})
                    reactions.append({
                        "id": rid_full,
                        "equation": detail.get('EQUATION'),
                        "enzymes": detail.get('ENZYME'),
                        "pathways": detail.get('PATHWAY'),
                    })
            except Exception as e:
                reactions = [{"id": rid, "equation": f"Error: {e}"} for rid in reaction_ids]

        ec_res = self.enzymes_of_compound(compound_id)
        enzymes = [it.get('id') for it in ec_res.get('data', [])][:max_items]

        # 3) 从反应映射 KO → 基因（可选物种）
        kos: List[str] = []
        for r in reaction_ids[:max_items]:
            k = self.kos_of_reaction(r)
            for it in k.get('data', [])[:max_items]:
                kos.append(it.get('id'))
        # 去重
        kos = list(dict.fromkeys(kos))[:max_items]

        genes: List[Dict[str, str]] = []
        for ko in kos:
            g = self.genes_of_ko(ko, organism=organism)
            genes.extend(g.get('data', [])[:max_items])
        # 基因去重
        seen = set(); genes_unique = []
        for it in genes:
            gid = it.get('id')
            if gid and gid not in seen:
                seen.add(gid); genes_unique.append(it)

        result = {
            "status": "success",
            "compound": {"id": f"cpd:{compound_id}", "label": resolved_label},
            "pathways": pathways,
            "reactions": reactions,
            "enzymes": enzymes,
            "kos": kos,
            "genes": genes_unique,
        }

        if return_envelope:
            entities: List[Dict[str, Any]] = []
            entities.append({"type": "compound", "id": f"cpd:{compound_id}", "label": resolved_label, "attrs": {}})
            entities.extend({"type": "pathway", "id": pid, "label": None, "attrs": {}} for pid in pathways)
            entities.extend({"type": "reaction", "id": r.get('id'), "label": None, "attrs": {"equation": r.get('equation')}} for r in reactions)
            entities.extend({"type": "enzyme", "id": eid, "label": None, "attrs": {}} for eid in enzymes)
            entities.extend({"type": "ko", "id": kid, "label": None, "attrs": {}} for kid in kos)
            entities.extend({"type": "gene", "id": g.get('id'), "label": g.get('value'), "attrs": {}} for g in genes_unique)
            return {
                "version": "bioflow-1.0",
                "op": "kegg.workflow.compound_to_gene",
                "query": {"name_or_id": name_or_id, "organism": organism, "max_items": max_items},
                "status": "success",
                "entities": list(entities),
                "summary": {
                    "n_pathways": len(pathways),
                    "n_reactions": len(reactions),
                    "n_enzymes": len(enzymes),
                    "n_kos": len(kos),
                    "n_genes": len(genes_unique)
                }
            }
        return result

if __name__ == "__main__":
    try:
        print(f"[RUNNING FILE] {__file__}")
    except Exception:
        pass

    parser = argparse.ArgumentParser(description="KEGG Tool Workflow CLI")
    parser.add_argument("--workflow", type=str, required=True, help="一键链路：输入化合物名或C编号，如 glucose 或 C00031")
    parser.add_argument("--org", type=str, default=None, help="限定物种（如 eco/hsa）")
    parser.add_argument("--limit", type=int, default=20, help="每步最大条目数")
    # The --envelope argument is no longer used for tool instantiation
    parser.add_argument("--envelope", action="store_true", help="此参数当前未使用")
    args = parser.parse_args()

    # [结构性修复] 现在直接实例化即可，无需传递参数
    tool = KeggTool()
    res = tool.resolve_compound_workflow(name_or_id=args.workflow,
                                         organism=args.org,
                                         max_items=args.limit,
                                         return_envelope=bool(args.envelope))
                                         
    print(json.dumps(res, ensure_ascii=False, indent=2))