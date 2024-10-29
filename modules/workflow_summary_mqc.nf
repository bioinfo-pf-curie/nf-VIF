process WORKFLOW_SUMMARY_MQC {
  input:
  val mqc_yaml_content

  output:
  path('workflow_summary_mqc.yaml')

  """
  echo -e "${mqc_yaml_content}" > workflow_summary_mqc.yaml
  """
}