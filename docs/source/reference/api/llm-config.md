---
myst:
  html_meta:
    description: Reference for the fields of an LLM configuration, the dicts that CyteType.run() accepts in llm_configs to use your own model provider.
---

# LLM configuration

`CyteType.run()` accepts your own models in `llm_configs`, a list of dicts. Each dict is
validated against `cytetype.api.LLMModelConfig` before any file is uploaded: a missing
required field raises a `pydantic.ValidationError`, and keys that are not fields of the
model are ignored, so check their spelling (field names are case-sensitive). The server
then checks every model before the job starts, and `run()` raises `LLMValidationError` if
one fails. [Use your own LLM](../guides/llm-configuration) explains when to use these fields,
with examples.

```python
adata = annotator.run(
    study_context="Human PBMCs from a healthy donor.",
    llm_configs=[
        {"provider": "openai", "name": "<model-name>", "apiKey": "<your-api-key>"},
    ],
)
```

Each configuration needs either `apiKey` or all three AWS fields (`awsAccessKeyId`,
`awsSecretAccessKey` and `awsDefaultRegion`).

```{eval-rst}
.. py:class:: cytetype.api.LLMModelConfig

   The fields of one entry in ``llm_configs``.

   .. py:attribute:: provider
      :type: str

      Required. The model provider: ``anthropic``, ``bedrock``, ``deepseek``,
      ``fireworks``, ``google``, ``groq``, ``huggingface``, ``mistral``, ``nebius``,
      ``openai``, ``openrouter``, ``together``, ``vertex`` or ``xai``.

   .. py:attribute:: name
      :type: str

      Required. The model name, as the provider names it.

   .. py:attribute:: apiKey
      :type: str | None
      :value: None

      Your API key for the provider. Required unless the three AWS fields are set.

   .. py:attribute:: baseUrl
      :type: str | None
      :value: None

      The base URL of the provider's API. The ``openai`` provider uses it to reach
      OpenAI-compatible servers, and ``vertex`` uses it too; other providers ignore it.

   .. py:attribute:: awsAccessKeyId
      :type: str | None
      :value: None

      AWS access key ID for ``bedrock``. Set all three AWS fields together.

   .. py:attribute:: awsSecretAccessKey
      :type: str | None
      :value: None

      AWS secret access key for ``bedrock``.

   .. py:attribute:: awsDefaultRegion
      :type: str | None
      :value: None

      AWS region for ``bedrock``.

   .. py:attribute:: modelSettings
      :type: dict[str, Any] | None
      :value: None

      Extra settings for the requests to the model, such as ``temperature``.

   .. py:attribute:: targetAgents
      :type: list[str] | None
      :value: []

      The CyteType agents that use this model: ``contextualizer``, ``annotator`` or
      ``reviewer``. An empty list means all agents. ``summarizer``, ``clinician`` and
      ``chat`` are also accepted but have no effect: the current workflow does not use
      the summarizer, and target analysis and Chat always use CyteType's hosted models.

   .. py:attribute:: skipValidation
      :type: bool
      :value: False

      Has no effect: the server validates every model regardless.

   .. py:attribute:: allowFallback
      :type: bool
      :value: False

      If true, an agent tries this model only after its models without
      ``allowFallback``.
```
